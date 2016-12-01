#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
CA_PATH=$1
ASM_DIR=$2;
ASM_PREFIX=$3;
NUM_THREADS=$4;
OVL_MER=$5;
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$CA_PATH:$MYPATH:$PATH;

set -e

if [ ! -e "$ASM_DIR/self_map.success" ];then
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d consensus > $ASM_DIR/unitigs.initial.fa && \
nucmer -l 41 -c 400 -b 400 -p $ASM_DIR/asm_to_asm  $ASM_DIR/unitigs.initial.fa  $ASM_DIR/unitigs.initial.fa && \
rm $ASM_DIR/unitigs.initial.fa && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/asm_to_asm.delta| delta-filter -q -o 20 /dev/stdin|show-coords -lcHr /dev/stdin | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl 10000 | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' > $ASM_DIR/duplicates.txt && \
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d layout | awk '{if($1 ~/^unitig/){unitig=$2;}else if($1~/^FRG/){print $5" utg"unitig}}' | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.txt'");while($l=<FILE>){chomp($l);$d{$l}=1}}{print $F[0],"\n" if(defined($d{$F[1]}));}' > $ASM_DIR/duplicates.iid.txt && \
touch $ASM_DIR/self_map.success
fi

if [ ! -e "$ASM_DIR/unitig_mer.success" ];then
tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 5 -U -d layout |tr -d '-' | \
awk 'BEGIN{print ">unique unitigs"}{if($1 == "cns"){seq=$2}else if($1 == "data.unitig_coverage_stat" && $2>=5){print seq"N"}}' | \
jellyfish count -L 10 -C -m $OVL_MER -s $ESTIMATED_GENOME_SIZE -t $NUM_THREADS -o $ASM_DIR/unitig_mers /dev/fd/0 && \
touch $ASM_DIR/unitig_mer.success
fi

if [ -e "$ASM_DIR/self_map.success" ] && [ -e "$ASM_DIR/unitig_mer.success" ];then
rm -rf $ASM_DIR/ovlStoreBackup
mkdir -p $ASM_DIR/ovlStoreBackup && \
mv $ASM_DIR/{4-unitigger,5-consensus,5-consensus-coverage-stat,5-consensus-insert-sizes,genome.tigStore,genome.ovlStore} $ASM_DIR/ovlStoreBackup
else
echo "Failed to self-map or build unitig mer database: no $ASM_DIR/self_map.success or $ASM_DIR/unitig_mer.success"
exit
fi

if [ ! -e "$ASM_DIR/overlap_filter.success" ];then
overlapStore -d $ASM_DIR/ovlStoreBackup/$ASM_PREFIX.ovlStore | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.iid.txt'");while($l=<FILE>){chomp($l);$diid{$l}=1}}{if($F[0]<$F[1]){if(not(defined($diid{$F[0]})) && not(defined($diid{$F[1]}))){ print join(" ",@F[0..6]),"\n"}}}'  | filter_overlap_file -t $NUM_THREADS <(gatekeeper  -dumpfragments -withsequence $ASM_DIR/$ASM_PREFIX.gkpStore| grep -P '^fragmentIdent|^fragmentSequence' | \
perl -ane 'BEGIN{$flag=1}{if($flag){print ">";}print "$F[2]\n";$flag=1-$flag;}') $ASM_DIR/unitig_mers /dev/fd/0 | convertOverlap -ovl | gzip -c > $ASM_DIR/overlaps_dedup.ovb.gz && touch $ASM_DIR/overlap_filter.success
fi

if [ ! -e "$ASM_DIR/overlapStore_rebuild.success" ];then
overlapStoreBuild -o $ASM_DIR/$ASM_PREFIX.ovlStore -M 65536 -g $ASM_DIR/$ASM_PREFIX.gkpStore $ASM_DIR/overlaps_dedup.ovb.gz 1>$ASM_DIR/overlapStore.rebuild.err 2>&1 && rm $ASM_DIR/overlaps_dedup.ovb.gz && touch $ASM_DIR/overlapStore_rebuild.success
fi

#overlapStore -d $ASM_DIR/ovlStoreBackup/$ASM_PREFIX.ovlStore | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.iid.txt'");while($l=<FILE>){chomp($l);$diid{$l}=1}}{if($F[0]<$F[1]){if(not(defined($diid{$F[0]})) && not(defined($diid{$F[1]}))){ print join(" ",@F[0..6]),"\n"}}}' | convertOverlap -ovl | gzip -c > $ASM_DIR/overlaps_dedup.ovb.gz && \
#overlapStoreBuild -o $ASM_DIR/$ASM_PREFIX.ovlStore -M 65536 -g $ASM_DIR/$ASM_PREFIX.gkpStore $ASM_DIR/overlaps_dedup.ovb.gz 1>$ASM_DIR/overlapStore.rebuild.err 2>&1

