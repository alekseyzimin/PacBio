#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
CA_PATH=$1
ASM_DIR=$2;
ASM_PREFIX=$3;
NUM_THREADS=$4;
OVL_MER=$5;
PLOIDY=$6;

if [ $PLOIDY -gt 1 ];then
MERGE_LEN=20000
HAP_SIM_RATE=90
else
MERGE_LEN=10000
HAP_SIM_RATE=95
fi

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$CA_PATH:$MYPATH:$PATH;

set -e

#-CS option in bogart is buggy -- I think it is still worth to place contains into singletons to turn them into regular unitigs for subsequent merging.  and then afterwards we can eliminate overlaps for those reads that end up in a unitig with only one parent
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 2 -U -d layout | perl -ane '{if($F[0] eq "unitig"){$utg="utg".$F[1];}elsif($F[0] eq "FRG" && $F[6]==0){$maximal{$utg}++;}}END{foreach $u(keys %maximal){print $u,"\n" if($maximal{$u}==1);}}' > $ASM_DIR/duplicates.txt

#here we map the unitigs against themselves to figure out which ones are redundant , and then record the reads in the redundant unitigs to eliminate their overlaps
if [ ! -e "$ASM_DIR/self_map.success" ];then
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d consensus |  awk '{if($1 ~/^>/){h=$1}else{if(length($1)>500) print length($1)" "h" "$1}}' | sort -S 50% -nrk1 | awk '{print $2"\n"$3}' | ufasta extract -v -f $ASM_DIR/duplicates.txt /dev/stdin >  $ASM_DIR/unitigs.ref.fa && \
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d consensus | ufasta extract -v -f $ASM_DIR/duplicates.txt /dev/stdin > $ASM_DIR/unitigs.qry.fa && \
nucmer -t $NUM_THREADS --batch $(($(stat -c%s "$ASM_DIR/unitigs.ref.fa")/25)) -l 31 -c 200 -b 100 -p $ASM_DIR/asm_to_asm  $ASM_DIR/unitigs.ref.fa  $ASM_DIR/unitigs.qry.fa && \
touch $ASM_DIR/self_map.success
fi

if [ ! -e "$ASM_DIR/filter_map.success" ] && [ -e $ASM_DIR/self_map.success ];then
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/asm_to_asm.delta > $ASM_DIR/asm_to_asm.noself.delta &&  \
parallel_delta-filter.sh $ASM_DIR/asm_to_asm.noself -q $NUM_THREADS && \
show-coords -lcHr -I $HAP_SIM_RATE $ASM_DIR/asm_to_asm.noself.fdelta | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl $MERGE_LEN | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' >> $ASM_DIR/duplicates.txt && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/asm_to_asm.delta| show-coords -lcH /dev/stdin | awk '{if($12>$13 && $10>int("'$HAP_SIM_RATE'") && $16>90) print $NF}' >> $ASM_DIR/duplicates.txt && \
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d layout | awk '{if($1 ~/^unitig/){unitig=$2;}else if($1~/^FRG/){print $5" utg"unitig}}' | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.txt'");while($l=<FILE>){chomp($l);$d{$l}=1}}{print $F[0],"\n" if(defined($d{$F[1]}));}' > $ASM_DIR/duplicates.iid.txt && \
rm -f $ASM_DIR/unitigs.{ref,qry}.fa $ASM_DIR/asm_to_asm.noself.{f,}delta && \
touch $ASM_DIR/filter_map.success
fi

#here we examine 22-mers in the unique unitigs and build a database of them that have a count of 9 or above
if [ ! -e "$ASM_DIR/unitig_mer.success" ];then
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d layout |tr -d '-' | \
awk 'BEGIN{print ">unique unitigs"}{if($1 == "cns"){seq=$2}else if($1 == "data.unitig_coverage_stat" && $2>=5){print seq"N"}}' | \
jellyfish count -L 7 -C -m $OVL_MER -s `tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d sizes | grep '^tigLenAssembled sum' | awk '{print $3}'` -t $NUM_THREADS -o $ASM_DIR/unitig_mers /dev/fd/0 && \
touch $ASM_DIR/unitig_mer.success
fi

#we backup the assembly
if [ -e "$ASM_DIR/filter_map.success" ] && [ -e "$ASM_DIR/unitig_mer.success" ];then
rm -rf $ASM_DIR/ovlStoreBackup
mkdir -p $ASM_DIR/ovlStoreBackup && \
mv $ASM_DIR/{4-unitigger,5-consensus,5-consensus-coverage-stat,5-consensus-insert-sizes,genome.tigStore,genome.ovlStore} $ASM_DIR/ovlStoreBackup
else
echo "Failed to self-map or build unitig mer database: no $ASM_DIR/self_map.success or $ASM_DIR/unitig_mer.success"
exit
fi

#we filter the overlaps -- if ALL k-mers in the overlap region are repetitive from above -- break overlap
if [ ! -e "$ASM_DIR/overlap_filter.success" ];then
overlapStore -d $ASM_DIR/ovlStoreBackup/$ASM_PREFIX.ovlStore | awk '{if($1<$2) print $0}' | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.iid.txt'");while($l=<FILE>){chomp($l);$diid[$l]=1}}{print join(" ",@F[0..6]),"\n" unless(defined($diid[$F[0]]) || defined($diid[$F[1]]));}'  | filter_overlap_file -t $NUM_THREADS <(gatekeeper  -dumpfragments -withsequence $ASM_DIR/$ASM_PREFIX.gkpStore| grep -P '^fragmentIdent|^fragmentSequence' | \
perl -ane 'BEGIN{$flag=1}{if($flag){print ">";}print "$F[2]\n";$flag=1-$flag;}') $ASM_DIR/unitig_mers /dev/fd/0 | convertOverlap -ovl | gzip -c > $ASM_DIR/overlaps_dedup.ovb.gz && touch $ASM_DIR/overlap_filter.success
fi

#rebuild the overlapstore
if [ ! -e "$ASM_DIR/overlapStore_rebuild.success" ];then
overlapStoreBuild -o $ASM_DIR/$ASM_PREFIX.ovlStore -M 65536 -g $ASM_DIR/$ASM_PREFIX.gkpStore $ASM_DIR/overlaps_dedup.ovb.gz 1>$ASM_DIR/overlapStore.rebuild.err 2>&1 && rm $ASM_DIR/overlaps_dedup.ovb.gz && touch $ASM_DIR/overlapStore_rebuild.success
fi

#overlapStore -d $ASM_DIR/ovlStoreBackup/$ASM_PREFIX.ovlStore | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.iid.txt'");while($l=<FILE>){chomp($l);$diid{$l}=1}}{if($F[0]<$F[1]){if(not(defined($diid{$F[0]})) && not(defined($diid{$F[1]}))){ print join(" ",@F[0..6]),"\n"}}}' | convertOverlap -ovl | gzip -c > $ASM_DIR/overlaps_dedup.ovb.gz && \
#overlapStoreBuild -o $ASM_DIR/$ASM_PREFIX.ovlStore -M 65536 -g $ASM_DIR/$ASM_PREFIX.gkpStore $ASM_DIR/overlaps_dedup.ovb.gz 1>$ASM_DIR/overlapStore.rebuild.err 2>&1

