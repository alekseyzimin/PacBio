#!/bin/bash
# MaSuRCA pipeline
# Copyright (C) 2017 by Alexey Zimin.
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#this is a script to run fragScaff, it should bw placed in the same folder as fragScaff.pl
#assumes that NCBI blastn, and bedtools are on the PATH

PATH=/software/Linux64/bedtools2/bin:/home/alekseyz/ncbi-blast-2.5.0+/bin:/genome2/raid/alekseyz/MaSuRCA/build/inst/bin:$PATH

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
#these parameters are optimozed for chromium defaults are 10000 and 20000
E=20000
o=50000
PREFIX="scaff"
let COUNT=100
NUM_THREADS=1
ASSEMBLY="assembly.fa"
FQ="barcoded.fastq"
j=1.25
u=2.0
DRYRUN=0

while [[ $# > 0 ]]
do
key="$1"
case $key in
-p|--prefix)
PREFIX="$2"
shift
;;
-t|--threads)
NUM_THREADS="$2"
shift
;;
-j)
j="$2"
shift
;;
-u)
u="$2"
shift
;;
-E)
E="$2"
shift
;;
-o)
o="$2"
shift
;;
-d|--dryrun)
DRYRUN="$2"
shift
;;
-c|--count)
COUNT="$2"
shift
;;
-a|--assembly)
ASSEMBLY="$2"
shift
;;
-b|--barcoded-fastq)
FQ="$2"
shift
;;
-v|--verbose)
set -x
;;
-h|--help|-u|--usage)
echo "Usage: fragScaff.sh -t <NUM_THREADS> -p <PREFIX default scaff> -a <ASSEMBLED_SCAFFOLDS.FA> -c <MIN_COUNT default 100> -b <BARCODED.FASTQ with path>"
exit 0
;;
*)
echo "Unknown option $1"
exit 1        # unknown option
;;
esac
shift
done

if [ ! -e $ASSEMBLY ];then
echo "ERROR assembly $ASSEMBLY does not exist"
echo "Usage: fragScaff.sh -t <NUM_THREADS> -p <PREFIX default scaff> -a <ASSEMBLED_SCAFFOLDS.FA> -c <MIN_COUNT default 100> -b <BARCODED.FASTQ with path>"
exit
fi

if [ ! -e $FQ ];then
echo "ERROR input barcoded.fastq $FQ does not exist"
echo "Usage: fragScaff.sh -t <NUM_THREADS> -p <PREFIX default scaff> -a <ASSEMBLED_SCAFFOLDS.FA> -c <MIN_COUNT default 100> -b <BARCODED.FASTQ with path>"
exit
fi

#prepare the assembly

if [ ! -e $PREFIX.prepare.success ];then
rm -f $PREFIX.blast.success 
rm -f $PREFIX.bwaindex.success 
echo -n "Preparing assembly " && date
perl -ane '{if(substr($F[0],0,1) eq ">"){if(length($seq)>=3000){print ">$rn\n$seq\n"}elsif(length($seq)>0){print STDERR ">$rn\n$seq\n"} $rn=substr($F[0],1);$seq=""}else{$seq.=$F[0]}}END{if(length($seq)>=3000){print ">$rn\n$seq\n"}elsif(length($seq)>0){print STDERR ">$rn\n$seq\n"}}' $ASSEMBLY 1>$PREFIX.assembly.ge3000.fa 2>$PREFIX.assembly.lt3000.fa && \
touch $PREFIX.prepare.success || exit 
fi

if [ ! -e $PREFIX.blast.success ];then
rm -f $PREFIX.bed.success
echo -n "Running blastn " && date
makeblastdb -in $PREFIX.assembly.ge3000.fa -dbtype nucl -out $PREFIX.db
blastn -word_size 36 -perc_identity 95 -outfmt 6 -db $PREFIX.db -query $PREFIX.assembly.ge3000.fa  -out $PREFIX.blastResult -num_threads $NUM_THREADS && touch $PREFIX.blast.success || exit
fi 

if [ ! -e $PREFIX.bed.success ];then
rm -f $PREFIX.bamParse.success
echo -n "Making bed files " && date
perl $MYPATH/blast_self_alignment_filter.pl $PREFIX.blastResult 95 > $PREFIX.blastResult.bed && \
sortBed -i $PREFIX.blastResult.bed > $PREFIX.blastResult.sorted.bed && \
mergeBed -i  $PREFIX.blastResult.sorted.bed > $PREFIX.repeats.bed && \
perl $MYPATH/fasta_make_Nbase_bed.pl $ASSEMBLY > $PREFIX.Nbase.bed && \
touch $PREFIX.bed.success || exit
fi

#index
if [ ! -e $PREFIX.bwaindex.success ];then
rm -f $PREFIX.map.success
echo -n "Running bwa index " && date
bwa index $PREFIX.assembly.ge3000.fa -p $PREFIX.bwa && touch $PREFIX.bwaindex.success || exit
fi

#create file with "count barcode"
if [ ! -e $PREFIX.count.success ];then
rm -f $PREFIX.map.success
echo -n "Filtering barcoded reads " && date
zcat -f $FQ | grep '^@' | perl -ane '{$count{$F[1]}++ if(defined($F[1]))}END{foreach $k(keys %count){print "$count{$k} $k\n"}}'  > $PREFIX.countBarcode.txt && \
perl -ane '{$h{$F[1]}=$F[0]}
END{open(FILE,"zcat -f '$FQ' |");while($line=<FILE>){chomp($line);($h,$b)=split(/\s+/,$line);$rn=substr($h,1);if($h{$b}>int("'$COUNT'")){my @lines=();my @b=split(":",$b);push(@lines,"@".$b[2].":".$rn."\n");for($i=0;$i<3;$i++){$line=<FILE>;push(@lines,$line);} print @lines if(length($lines[1])>=100);}else{for($i=0;$i<3;$i++){$line=<FILE>;}}}}' $PREFIX.countBarcode.txt > $PREFIX.readsWithBarcode.ge100.fastq && \
touch $PREFIX.count.success ||exit
fi

#map the reads
if [ ! -e $PREFIX.map.success ];then
rm -f $PREFIX.bamParse.success 
rm -f $PREFIX.reheader.success
echo -n "Mapping barcoded reads " && date
bwa mem -p -t 64 $PREFIX.bwa $PREFIX.readsWithBarcode.ge100.fastq 2>bwasterr | samtools view -bhS /dev/stdin | samtools sort -m64000000000 /dev/stdin $PREFIX.alignSorted && \
touch $PREFIX.map.success ||exit
fi

#reheader
if [ ! -e $PREFIX.reheader.success ];then
rm -f $PREFIX.bamParse.success
echo -n "Reheader bam file " && date
samtools view -H $PREFIX.alignSorted.bam > $PREFIX.header.sam && \
awk '{if($1>=100) {split($2,a,":");print "@RG\tID:"a[3]}}' $PREFIX.countBarcode.txt >> $PREFIX.header.sam && \
samtools reheader $PREFIX.header.sam $PREFIX.alignSorted.bam > $PREFIX.alignSorted.reheader.bam && \
touch $PREFIX.reheader.success ||exit
fi

#run fragScaff
#parse BAM file
if [ ! -e $PREFIX.bamParse.success ];then
rm -f $PREFIX.links.success
rm -f $PREFIX.alignSorted.reheader.bam.E20000.o50000.J.N.bamParse.log
echo -n "Parsing bam file " && date
$MYPATH/fragScaff.pl -B $PREFIX.alignSorted.reheader.bam -b 1 -J $PREFIX.repeats.bed -N $PREFIX.Nbase.bed -E $E -o $o && \
touch $PREFIX.bamParse.success ||exit
fi 

if [ ! -e $PREFIX.links.success ];then
rm -f $PREFIX.scaffold.success
rm -f $PREFIX.alignSorted.reheader.bam.E20000.o50000.J.link.fragScaff.log
echo -n "Creating links " && date
$MYPATH/fragScaff.pl -B $PREFIX.alignSorted.reheader.bam.E${E}.o${o}.J.N.bamParse -A -O $PREFIX.alignSorted.reheader.bam.E${E}.o${o}.J.link -t $NUM_THREADS && \
touch $PREFIX.links.success ||exit
fi

if [ ! -e $PREFIX.scaffold.success ];then
echo -n "Computing scaffolds " && date
rm -f $PREFIX.fragScaff.log
if [ $DRYRUN -le 0 ];then
$MYPATH/fragScaff.pl -B $PREFIX.alignSorted.reheader.bam.E${E}.o${o}.J.N.bamParse -K $PREFIX.alignSorted.reheader.bam.E${E}.o${o}.J.N.r1.links.txt -j $j -u $u -F $PREFIX.assembly.ge3000.fa -O $PREFIX && \
cat $PREFIX.assembly.lt3000.fa >> $PREFIX.fragScaff.assembly.fasta && \
echo "Scaffolding success, output in $PREFIX.fragScaff.assembly.fasta" && cat $PREFIX.fragScaff.log && \
touch $PREFIX.scaffold.success ||exit
else
echo "Dry run to optimize parameters -- no scaffolds output!!!"
$MYPATH/fragScaff.pl -B $PREFIX.alignSorted.reheader.bam.E${E}.o${o}.J.N.bamParse -K $PREFIX.alignSorted.reheader.bam.E${E}.o${o}.J.N.r1.links.txt -j $j -u $u -I -O $PREFIX && cat $PREFIX.fragScaff.log 
fi
fi

