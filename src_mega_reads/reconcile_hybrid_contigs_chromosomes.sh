#!/bin/bash
#this code aims at reconciling the hybrid contigs and the chromosomes of the previously produces assembly
#arguments are: reference chromosomes, hybrid contigs, hybrid posmap (frgctg), filtered delta-file of alignments of ref to hyb
#MUST HAVE MaSURCA bin on the PATH
REF_CHR=$1
HYB_CTG=$2
HYB_POS=$3
REF_HYB_FDELTA=$4

#function cleanup {
#if [ $$ -eq $pid ];then
#rm -rf .rerun 
#fi
#}
#trap cleanup EXIT SIGHUP SIGINT SIGTERM


set -x
#compute coverage from the posmap file
if [ ! -s $HYB_POS.coverage ] || [ -e .rerun ];then 
awk '{print $1" "$2" "$3"\n"$1" "$2" "$4}' $HYB_POS | grep -v F |grep -v R | sort -nk2 -k3n -S 10% | compute_coverage.pl > $HYB_POS.coverage
touch .rerun
fi

if [ ! -s  gap_coordinates.txt ] || [ -e .rerun ];then
splitFileAtNs $REF_CHR 1 > $REF_CHR.split.fasta
perl -ane '{$h{substr($F[1],3)}=$F[0]}END{while($line=<STDIN>){chomp($line);@f=split(/\s+/,$line);print "$f[0] $h{$f[1]} ",$f[2]+1," $f[3] $f[4]\n";}}' scaffNameTranslations.txt < genome.posmap.ctgscf | awk 'BEGIN{pg=0}{print $2" "pg" "$3;pg=$4}' > gap_coordinates.txt
touch .rerun
fi

if [ ! -s $REF_HYB_FDELTA.coords ] || [ -e .rerun ];then
show-coords -lcHr $REF_HYB_FDELTA | \
merge_matches_and_tile_coords_file.pl 250000 | \
merge_matches_and_tile_coords_file.pl 100000 | \
awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
awk '{if($16>5 || $7>5000 ) print $0}' > $REF_HYB_FDELTA.coords
touch .rerun
fi

#find split contigs
awk '{if($4<$5) print $4" "$5" "($4+$5)/2" "$NF" "$13; else print $5" "$4" "($4+$5)/2" "$NF" "$13;}' $REF_HYB_FDELTA.coords| \
sort -k4 -k1n | \
uniq -D -f 3 | \
awk '{if($NF != prev){offset=$2;prev=$NF;print $0}else if($2>offset){print $0;offset=$2;}else{print "contained "$0}}' > $REF_HYB_FDELTA.coords.split_contain

#create breaks files
grep -v ^contained $REF_HYB_FDELTA.coords.split_contain | uniq -D -f 3 | awk '{if($4==ctg){if($1>5000 && $1<$NF-5000) print "alnbreak "substr($4,4)" "$1" 0"}else{ctg=$4;if($2>5000 && $2<$NF-5000) print "alnbreak "substr($4,4)" "$2" 0"}}' > $REF_HYB_FDELTA.breaks

cat <(perl -ane 'chop;chop;print;print "\n"' A.tau\ Potential\ chimeric\ point.txt |awk '{print "break "substr($1,4)" "$2" "$3}') asm.posmap.ctgscf|grep -v ID | sort -nk2 -k3n -S 10% | awk '{if($1~ /^break/){if(dir == "f"){if($3-soffset >1000 && $3-soffset<len-1000)print "break_scf"$2"_"$3"_"$4" "ctg" "scf" "$3-soffset" 0";}else{if(eoffset-$3 >1000 && eoffset-$3<len-1000) print "break_scf"$2"_"$3"_"$4" "ctg" "scf" "eoffset-$3" 0";}}else{ctg=$1;soffset=$3;eoffset=$4; len=$4-$3;dir=$5;}}' > nanomap_breaks.txt

if [ ! -s $HYB_POS.coverage.w_breaks ] || [ -e .rerun ];then
cat nanomap_breaks.txt $REF_HYB_FDELTA.breaks $HYB_POS.coverage  | sort -nk2 -k3n -S 10% > $HYB_POS.coverage.w_breaks
touch .rerun
fi

grep -C 150 break $HYB_POS.coverage.w_breaks  | evaluate_splits.pl <(ufasta sizes -H $HYB_CTG | tr -d 'ctg') | sort -nk3 -S 10% >  $HYB_POS.coverage.w_breaks.validated

if [ ! -s $HYB_CTG.broken ] || [ -e .rerun ];then
break_contigs.pl <(grep -v "end" $HYB_POS.coverage.w_breaks.validated |awk '{if($4<4) print $0}') < $HYB_CTG > $HYB_CTG.broken
fi

#now we re-align the broken contigs to the reference
if [ ! -s $REF_CHR.$HYB_CTG.broken.delta ] || [ -e .rerun ];then
nucmer -p $REF_CHR.$HYB_CTG.broken -c 250 -l 31 $REF_CHR $HYB_CTG.broken
fi
if [ ! -s $REF_CHR.$HYB_CTG.broken.1.delta ] || [ -e .rerun ];then
delta-filter -1 -o 20 -i 98 $REF_CHR.$HYB_CTG.broken.delta > $REF_CHR.$HYB_CTG.broken.1.delta
fi

#now we merge/rebuild chromosomes
show-coords -lcHr $REF_CHR.$HYB_CTG.broken.1.delta | \
merge_matches_and_tile_coords_file.pl 250000 | \
merge_matches_and_tile_coords_file.pl 100000 | \
awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
awk '{if($16>5 || $7>5000 ) print $0}' > $REF_CHR.$HYB_CTG.broken.1.coords

# here we split everything so the contigs are "perfect"
cat $REF_CHR.$HYB_CTG.broken.1.coords  |  extract_single_best_match_coords_file.pl |grep Chr |reconcile_matches.pl gap_coordinates.txt  > reconciled_coords.txt

cat reconciled_coords.txt  | output_reconciled_scaffolds.pl $HYB_CTG.broken > $REF_CHR.$HYB_CTG.reconciled.fa

#grep ' break' validated_breaks.txt |perl -ane '{for($i=4;$i<=$#F;$i++){if($F[$i]=~/^break/){@f=split("_",$F[$i]);$report{"$f[1] $f[2]"}="contig $F[1] position $F[2] coverage $F[3]";last}}}END{open(FILE,"A.tau Potential chimeric point.txt");while($line=<FILE>){chop($line);chop($line);print $line;@f=split(/\s+/,$line);print " ",$report{"$f[0] $f[1]"},"\n";}}' > A.tau\ Potential\ chimeric\ point.txt.validated

