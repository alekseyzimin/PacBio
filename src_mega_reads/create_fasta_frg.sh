#!/bin/bash

COORDS=$1;
KMER=$2;
PB_SEQ=$3
GOOD_PB=$4
BAD_PB=$COORDS.bad_pb.txt

join_mega_reads_trim.pl $PB_SEQ ${COORDS}.allowed $KMER $BAD_PB $GOOD_PB < ${COORDS}.all.txt > $COORDS.fa;
fasta2frg.pl mr  < $COORDS.fa > $COORDS.frg;
/home/alekseyz/myprogs/getNumBasesPerReadInFastaFile.perl  $COORDS.fa   | awk '{n+=$1;m++}END{print n" "m" "n/m}'
