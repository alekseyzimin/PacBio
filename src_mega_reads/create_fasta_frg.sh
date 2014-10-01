#!/bin/bash
LEN=$1;
COORDS=$2;
KMER=$3;
PB_SEQ=$4
GOOD_PB=$5;

join_mega_reads_trim.pl $PB_SEQ $LEN ${COORDS}.allowed $KMER $GOOD_PB < ${COORDS}.all.txt > $COORDS.$LEN.fa;
fasta2frg.pl mr  < $COORDS.$LEN.fa > $COORDS.$LEN.frg;
/home/alekseyz/myprogs/getNumBasesPerReadInFastaFile.perl  $COORDS.$LEN.fa   | awk '{n+=$1;m++}END{print n" "m" "n/m}'
