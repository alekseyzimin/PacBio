#!/bin/bash
WORKDIR=$1  #input/output folder
FILENAME=$2  #filename of a coores file for matches of SINGLE pacbio read to the super-reads

#will fail if there is a gap in pacBio read coverage

./filter_matches.pl < $WORKDIR/$FILENAME  > $WORKDIR/$FILENAME.filtered
wc -l $WORKDIR/$FILENAME.filtered
cat $WORKDIR/$FILENAME.filtered | ./findPathAcrossPacbioRead.perl| dot  -T svg -o $WORKDIR/$FILENAME.filtered.svg
cat $WORKDIR/$FILENAME.filtered | ./build_mega_reads.pl 1>$WORKDIR/$FILENAME.debug.out 2>$WORKDIR/$FILENAME.megareads; 
grep ^used $WORKDIR/$FILENAME.debug.out | perl -ane '{print join(" ",@F[1..$#F]),"\n"}' > $WORKDIR/$FILENAME.filtered2
cat $WORKDIR/$FILENAME.filtered2 | ./findPathAcrossPacbioRead.perl| dot  -T svg -o $WORKDIR/$FILENAME.filtered2.svg
