#!/bin/bash

FILENAME=$1  #filename of a coords file for matches of SINGLE pacbio read to the super-reads
NAMESEQFILE=$2 #file containing pb read "name sequence"
EXEPATH=`dirname $0`

#will fail if there is a gap in pacBio read coverage

$EXEPATH/filter_matches.pl < $FILENAME  > $FILENAME.filtered
wc -l $FILENAME.filtered
#cat $FILENAME.filtered | $EXEPATH/findPathAcrossPacbioRead.perl| dot  -T svg -o $FILENAME.filtered.svg
cat $FILENAME.filtered | $EXEPATH/build_mega_reads.pl 1> $FILENAME.debug.out 2>$FILENAME.megareads;
#cat $FILENAME.filtered | $EXEPATH/build_mega_reads_kmerscore.pl 1> $FILENAME.debug.out 2>$FILENAME.megareads; 
grep ^used $FILENAME.debug.out | perl -ane '{print join(" ",@F[1..$#F]),"\n"}' > $FILENAME.filtered2
#cat $FILENAME.filtered2 | $EXEPATH/findPathAcrossPacbioRead.perl| dot  -T svg -o $FILENAME.filtered2.svg
if [ -L /genome7/raid/alekseyz/PB_ScerW303/assembly/work1 ];then
ln -s /genome7/raid/alekseyz/PB_ScerW303/assembly/work1
fi

if [ -s $FILENAME.megareads ];then
/home/alekseyz/myprogs/MaSuRCA/build/inst/bin/createFastaSuperReadSequences work1 <(awk '{print $1" "$5}' $FILENAME.megareads) -seqdiffmax 0 -min-ovl-len 69 -minreadsinsuperread 1  -good-sr-filename $FILENAME.megareads.names  -kunitigsfile /genome7/raid/alekseyz/PB_ScerW303/assembly/guillaumeKUnitigsAtLeast32bases_all.fasta -good-sequence-output-file $FILENAME.megareads.fa -super-read-name-and-lengths-file $FILENAME.megareads.sizes   2> work1/createFastaSuperReadSequences.errors.txt
fi
touch $FILENAME.megareads.fa
paste $FILENAME.megareads <(grep -v '>' $FILENAME.megareads.fa) > $FILENAME.megareads.wseq
sort -nrk4,4 $FILENAME.megareads.wseq | $EXEPATH/reconciliate_mega_reads.pl $NAMESEQFILE $EXEPATH | sort -nk2,2 > $FILENAME.megareads.sorted.wseq

