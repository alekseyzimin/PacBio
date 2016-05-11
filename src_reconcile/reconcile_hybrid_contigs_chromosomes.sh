#!/bin/bash
#this code aims at reconciling the hybrid contigs and the chromosomes of the previously produces assembly
#arguments are: reference chromosomes, hybrid contigs, hybrid posmap (frgctg), filtered delta-file of alignments of ref to hyb
#MUST HAVE MaSURCA bin on the PATH
REF_CHR=$1
CA_TERM_DIR=$2

set -xe

if [ ! -s  genome.scfdeg.fasta ] ||  [ ! -s  genome.posmap.frgscfdeg ];then
cat $CA_TERM_DIR/genome.{scf,deg}.fasta > genome.scfdeg.fasta
cat $CA_TERM_DIR/genome.posmap.frgscf.sorted $CA_TERM_DIR/genome.posmap.frgdeg  > genome.posmap.frgscfdeg
touch .rerun
fi

HYB_CTG="genome.scfdeg.fasta"
HYB_POS="genome.posmap.frgscfdeg"

#function cleanup {
#if [ $$ -eq $pid ];then
#rm -rf .rerun 
#fi
#}
#trap cleanup EXIT SIGHUP SIGINT SIGTERM

if [ ! -s $REF_CHR.$HYB_CTG.delta ];then
nucmer -p $REF_CHR.$HYB_CTG -c 250 -l 31 $REF_CHR $HYB_CTG
touch .rerun
fi

if [ ! -s $REF_CHR.$HYB_CTG.1.delta ] || [ -e .rerun ];then
delta-filter -1 -i 98 -o 20 $REF_CHR.$HYB_CTG.delta >$REF_CHR.$HYB_CTG.1.delta
touch .rerun
fi

#compute coverage from the posmap file
if [ ! -s $HYB_POS.coverage ];then 
awk '{print $1" "$2" "$3"\n"$1" "$2" "$4}' $HYB_POS | grep -v F |grep -v R | sort -nk2 -k3n -S 10% | compute_coverage.pl > $HYB_POS.coverage
touch .rerun
fi

if [ ! -s  gap_coordinates.txt ];then
splitFileAtNs $REF_CHR 1 > $REF_CHR.split.fasta
perl -ane '{$h{substr($F[1],3)}=$F[0]}END{while($line=<STDIN>){chomp($line);@f=split(/\s+/,$line);print "$f[0] $h{$f[1]} ",$f[2]+1," $f[3] $f[4]\n";}}' scaffNameTranslations.txt < genome.posmap.ctgscf | awk 'BEGIN{pg=0}{print $2" "pg" "$3;pg=$4}' > gap_coordinates.txt
touch .rerun
fi

if [ ! -s $REF_CHR.$HYB_CTG.1.coords ] || [ -e .rerun ];then
show-coords -lcHr $REF_CHR.$HYB_CTG.1.delta | \
merge_matches_and_tile_coords_file.pl 250000 | \
merge_matches_and_tile_coords_file.pl 100000 | \
awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
awk '{if($16>5 || $7>5000 ) print $0}' > $REF_CHR.$HYB_CTG.1.coords
touch .rerun
fi

#find split contigs
awk '{if($4<$5) print $4" "$5" "($4+$5)/2" "$NF" "$13; else print $5" "$4" "($4+$5)/2" "$NF" "$13;}' $REF_CHR.$HYB_CTG.1.coords| \
sort -k4 -k1n | \
uniq -D -f 3 | \
awk '{if($NF != prev){offset=$2;prev=$NF;print $0}else if($2>offset){print $0;offset=$2;}else{print "contained "$0}}' > $REF_CHR.$HYB_CTG.1.coords.split_contain

#create breaks files
grep -v ^contained $REF_CHR.$HYB_CTG.1.coords.split_contain | uniq -D -f 3 | awk '{if($4==ctg){if($1>5000 && $1<$NF-5000) print "alnbreak "substr($4,4)" "$1" 0"}else{ctg=$4;if($2>5000 && $2<$NF-5000) print "alnbreak "substr($4,4)" "$2" 0"}}' > $REF_CHR.$HYB_CTG.1.coords.breaks

if [ ! -s $HYB_POS.coverage.w_breaks ] || [ -e .rerun ];then
cat $REF_CHR.$HYB_CTG.1.coords.breaks $HYB_POS.coverage  | sort -nk2 -k3n -S 10% > $HYB_POS.coverage.w_breaks
touch .rerun
fi

grep -C 10 break $HYB_POS.coverage.w_breaks  | evaluate_splits.pl <(ufasta sizes -H $HYB_CTG | awk '{print substr($1,4)" "$2}') | sort -nk3 -S 10% >  $HYB_POS.coverage.w_breaks.validated

if [ ! -s $HYB_CTG.broken ] || [ -e .rerun ];then
break_contigs.pl <(grep -v "end" $HYB_POS.coverage.w_breaks.validated |awk '{if($4<2) print $0}') < $HYB_CTG > $HYB_CTG.broken
touch .rerun
fi

#now we re-align the broken contigs to the reference
if [ ! -s $REF_CHR.$HYB_CTG.broken.delta ] || [ -e .rerun ];then
nucmer -p $REF_CHR.$HYB_CTG.broken -c 250 -l 31 $REF_CHR $HYB_CTG.broken
touch .rerun
fi
if [ ! -s $REF_CHR.$HYB_CTG.broken.1.delta ] || [ -e .rerun ];then
delta-filter -1 -o 20 -i 98 $REF_CHR.$HYB_CTG.broken.delta > $REF_CHR.$HYB_CTG.broken.1.delta
touch .rerun
fi

rm -rf .rerun

#now we merge/rebuild chromosomes
show-coords -lcHr $REF_CHR.$HYB_CTG.broken.1.delta | \
merge_matches_and_tile_coords_file.pl 250000 | \
merge_matches_and_tile_coords_file.pl 100000 | \
awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
awk '{if($16>5 || $7>5000 ) print $0}' > $REF_CHR.$HYB_CTG.broken.1.coords

# here we split everything so the contigs are "perfect"
cat $REF_CHR.$HYB_CTG.broken.1.coords  |  extract_single_best_match_coords_file.pl  |reconcile_matches.pl gap_coordinates.txt  > reconciled_coords.txt

cat reconciled_coords.txt  | output_reconciled_scaffolds.pl $HYB_CTG.broken > $REF_CHR.$HYB_CTG.reconciled.fa

