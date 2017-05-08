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
parallel_delta-filter.sh $DELTAFILE '-1 -l 100 -o 10 ' 9 && mv $DELTAFILE.fdelta $DELTAFILE.1.delta && \
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

#here we map the contigs against themselves to figure out which ones are redundant
if [ ! -e self_map.contigs.success ];then
rm -f filter_map.contigs.success
perl -ane 'BEGIN{$seq=""}{if($F[0] =~ /^\>/){if(length($seq)>500){print length($seq)," $rn $seq\n";}$seq="";$rn=$F[0];}else{$seq.=$F[0];}}END{print length($seq)," $rn $seq\n";}' $REFN.$QRYN.all.polished.fa | sort -S 50% -nrk1 | awk '{print $2"\n"$3}' >  scaffolds.ref.fa && \
nucmer -t $NUM_THREADS --batch $(($(stat -c%s "scaffolds.ref.fa")/25)) -l 51 -c 100 -b 100 -p sasm_to_sasm  scaffolds.ref.fa  $REFN.$QRYN.all.polished.fa && \
touch self_map.contigs.success || exit
fi

if [ ! -e filter_map.contigs.success ];then
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' sasm_to_sasm.delta > sasm_to_sasm.noself.delta && \
delta-filter -i 98 -q -o 10 sasm_to_sasm.noself.delta > sasm_to_sasm.noself.fdelta && \
show-coords -lcHr sasm_to_sasm.noself.fdelta | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl 500 | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' > sduplicates.txt && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' sasm_to_sasm.delta| show-coords -lcH /dev/stdin | awk '{if($12>$13 && $10>98 && $16>90) print $NF}' >> sduplicates.txt && \
ufasta extract -v -f sduplicates.txt $REFN.$QRYN.all.polished.fa > $REFN.$QRYN.all.polished.deduplicated.fa && \
touch filter_map.contigs.success || exit
fi

echo "Output sequences in $REFN.$QRYN.all.polished.deduplicated.fa"
ufasta n50 -A -S -N50 $REFN.$QRYN.all.polished.deduplicated.fa
