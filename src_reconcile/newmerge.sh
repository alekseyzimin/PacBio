#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
REF=$1
QRY=$2
NUM_THREADS=$3

REFN=`basename $REF`
QRYN=`basename $QRY`
DELTAFILE=$REFN.$QRYN


#nucmer
if [ ! -e align.success ];then
nucmer -t $NUM_THREADS -p  $DELTAFILE -l 51 -c 200 $REF $QRY && touch align.success || exit
fi

#delta-filter
if [ ! -e filter.success ];then
rm -f show.success
rm -f add_not_alingning.success
rm -f merge.success
parallel_delta-filter.sh $DELTAFILE '-r -l 1000' 9 && mv $DELTAFILE.fdelta $DELTAFILE.r.delta && \
parallel_delta-filter.sh $DELTAFILE '-1 -l 100' 9 && mv $DELTAFILE.fdelta $DELTAFILE.1.delta && \
touch filter.success || exit
fi

if [ ! -e merge.success ];then
show-coords -lcHq -I 99 -L 3000 $DELTAFILE.r.delta | extract_merges.pl $QRY > merges.txt && merge_contigs.pl < merges.txt| create_merged_sequences.pl $REF merges.txt > $REFN.$QRYN.merged.fa && touch merge.success || exit
fi

if [ ! -e add_not_aligning.success ];then
#add the sequences that did not align
ufasta extract -v -f <(show-coords -lcH -I 99 $DELTAFILE.1.delta| perl -ane '{$palign{$F[-1]}+=$F[-4];}END{foreach $k(keys %palign){print $k,"\n" if($palign{$k}>20)}}') $QRY > $REFN.$QRYN.extra.fa && \
cat $REFN.$QRYN.merged.fa $REFN.$QRYN.extra.fa > $REFN.$QRYN.all.fa && \
touch add_not_aligning.success || exit
fi

if [ ! -e replace_consensus.success ];then
show-coords -lcHr $DELTAFILE.1.delta | reconcile_consensus.pl $REFN.$QRYN.all.fa $QRY > $REFN.$QRYN.all.polished.fa && \
touch replace_consensus.success ||exit
fi


echo "Output sequences in $REFN.$QRYN.all.polished.fa"
ufasta n50 -A -S -N50 $REFN.$QRYN.all.polished.fa
