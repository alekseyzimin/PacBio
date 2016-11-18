#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
CA_PATH=$1
ASM_DIR=$2;
ASM_PREFIX=$3;
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$CA_PATH:$MYPATH:$PATH;

set -e

tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 3 -U -d consensus > $ASM_DIR/unitigs.initial.fa

nucmer -l 51 -c 400 -b 400 -p $ASM_DIR/asm_to_asm  $ASM_DIR/unitigs.initial.fa  $ASM_DIR/unitigs.initial.fa

awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/asm_to_asm.delta| delta-filter -q -o 20 /dev/stdin|show-coords -lcHr /dev/stdin | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl 10000 | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' > $ASM_DIR/duplicates.txt

tigStore -g $ASM_DIR/$ASM_PREFIX.gkpStore -t $ASM_DIR/$ASM_PREFIX.tigStore 3 -U -d layout | awk '{if($1 ~/^unitig/){unitig=$2;}else if($1~/^FRG/){print $5" utg"unitig}}' | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.txt'");while($l=<FILE>){chomp($l);$d{$l}=1}}{print $F[0],"\n" if(defined($d{$F[1]}));}' > $ASM_DIR/duplicates.iid.txt

overlapStore -d $ASM_DIR/$ASM_PREFIX.ovlStore | perl -ane 'BEGIN{open(FILE,"'$ASM_DIR/duplicates.iid.txt'");while($l=<FILE>){chomp($l);$diid{$l}=1}}{if(not(defined($diid{$F[0]})) && not(defined($diid{$F[1]}))){ print join(" ",@F[0..6]),"\n"}}'  | convertOverlap -ovl | gzip -c > $ASM_DIR/overlaps.ovb.gz && \
mkdir -p $ASM_DIR/ovlStoreBackup && \
mv $ASM_DIR/{4-unitigger,5-consensus,5-consensus-coverage-stat,5-consensus-insert-sizes,genome.tigStore,genome.ovlStore} $ASM_DIR/ovlStoreBackup && \
overlapStoreBuild -o $ASM_DIR/$ASM_PREFIX.ovlStore -M 65536 -g $ASM_DIR/$ASM_PREFIX.gkpStore $ASM_DIR/overlaps.ovb.gz 1>$ASM_DIR/overlapStore.rebuild.err 2>&1

