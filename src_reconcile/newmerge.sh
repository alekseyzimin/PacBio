#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
REF=$1
QRY=$2
DELTAFILE=$3

REFN=`basename $REF`
QRYN=`basename $QRY`

#delta-filter
if [ ! -e filter.success ];then
rm -f show.success
rm -f add_not_alingning.success
rm -f merge.success
parallel_delta-filter.sh $DELTAFILE '-r -l 1000' 9 && mv $DELTAFILE.fdelta $DELTAFILE.r.delta && \
parallel_delta-filter.sh $DELTAFILE '-1 -l 1000' 9 && mv $DELTAFILE.fdelta $DELTAFILE.1.delta && \
touch filter.success || exit
fi

if [ ! -e merge.success ];then
show-coords -lcHq -I 99 -L 3000 $DELTAFILE.r.delta | extract_merges.pl $QRY > merges.txt && merge_contigs.pl < merges.txt| create_merged_sequences.pl $REF merges.txt > $REFN.$QRYN.merged.fa && touch merge.success || exit
fi

if [ ! -e add_not_aligning.success ];then
#add the sequences that did not align
ufasta extract -v -f <(show-coords -lcH -I 99 -L 1000 $DELTAFILE.1.delta| perl -ane '{$palign{$F[-1]}+=$F[-4];}END{foreach $k(keys %palign){print $k,"\n" if($palign{$k}>20)}}') $QRY > $REFN.$QRYN.extra.fa && \
cat $REFN.$QRYN.merged.fa $REFN.$QRYN.extra.fa > $REFN.$QRYN.all.fa && \
touch add_not_aligning.success || exit
fi

echo "Output sequences in $REFN.$QRYN.all.fa"
ufasta n50 -A -S -N50 $REFN.$QRYN.all.fa
