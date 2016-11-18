#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

ASM_DIR=$1;
ASM_PREFIX=$2;

nucmer -l 31 -c 400 -b 400 -p asm_to_asm $ASM_DIR/9-terminator/$ASM_PREFIX.ctg.fasta $ASM_DIR/9-terminator/$ASM_PREFIX.ctg.fasta
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' asm_to_asm.delta| delta-filter -q -o 20 /dev/stdin|show-coords -lcHr /dev/stdin | awk '{if($12>$13) print $0}' |/genome2/raid/alekseyz/MaSuRCA/build/inst/bin/merge_matches_and_tile_coords_file.pl 100000 | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' > $ASM_DIR/9-terminator/duplicates.txt
ufasta extract -f $ASM_DIR/9-terminator/duplicates.txt $ASM_DIR/9-terminator/$ASM_PREFIX.ctg.fasta > $ASM_DIR/9-terminator/$ASM_PREFIX.ctg.alt.fasta
ufasta extract -v -f $ASM_DIR/9-terminator/duplicates.txt $ASM_DIR/9-terminator/$ASM_PREFIX.ctg.fasta > $ASM_DIR/9-terminator/$ASM_PREFIX.ctg.main.fasta
exit


perl -ane '{$d{substr($F[0],3)}=1}END{open(FILE,"'$ASM_DIR/9-terminator/$ASM_PREFIX.posmap.frgctg'");while($l=<FILE>){chomp($l);@f=split(/\s+/,$l);print $f[0],"\n" if(defined($d{$f[1]}));}}' $ASM_DIR/9-terminator/duplicates.uid.txt
grep ^FRG $ASM_DIR/9-terminator/$ASM_PREFIX.iidtouid | perl -ane '{$iid{$F[2]}=$F[1]}END{open(FILE,"'$ASM_DIR/9-terminator/duplicates.uid.txt'");while($l=<FILE>){chomp($l);print $iid{$l},"\n";}}' > $ASM_DIR/9-terminator/duplicates.iid.txt
overlapStore -d $ASM_DIR/$ASM_PREFIX.ovlStore | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/9-terminator/duplicates.uid.txt'");while($l=<FILE>){chomp($l);$diid{$l}=1}}{if(not(defined($diid{$F[0]})) && not(defined($diid{$F[1]}))){ print join(" ",@F[0..6])}}'  | \
convertOverlap -ovl | gzip -c > $ASM_DIR/overlaps.ovb.gz &&  \
rm -rf $ASM_DIR/{4-unitigger,5-consensus,5-consensus-coverage-stat,5-consensus-insert-sizes,genome.tigStore} && \
mkdir -p $ASM_DIR/ovlStoreBackup && \
mv  $ASM_DIR/$ASM_PREFIX.ovlStore $ASM_DIR/ovlStoreBackup && \
overlapStoreBuild -o $ASM_DIR/$ASM_PREFIX.ovlStore -M 16384 -g $ASM_DIR/$ASM_PREFIX.gkpStore $ASM_DIR/overlaps.ovb.gz 1>$ASM_DIR/overlapStore.rebuild.err 2>&1

