#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
MERGE=0
SIMILARITY_RATE=97

function error_exit {
    echo "$1" >&2
    exit "${2:-1}"
}


#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            NUM_THREADS="$2"
            shift
            ;;
        -q|--query)
            QRY="$2"
            shift
            ;;
        -i|--identity)
            SIMILARITY_RATE="$2"
            shift
            ;;
        -m|--merge-slack)
            MERGE="$2"
            shift
            ;;
        -r|--reference)
            REF="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: polish_with_illumina_assembly.sh -r <sequence to be polished> -q <sequence to polish with> -t <number of threads> -i <minimum sequence identity percentage> -m <merge polishing sequence alignments slack (advanced)> "
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [ ! -s $REF ];then
error_exit "reference $REF does not exist or size zero"
fi

if [ ! -s $QRY ];then
error_exit "polishing sequence $QRY does not exist or size zero"
fi

REFN=`basename $REF`
QRYN=`basename $QRY`
DELTAFILE=$REFN.$QRYN

#nucmer
if [ ! -e polish_align.success ];then
nucmer -t $NUM_THREADS -p  $DELTAFILE -c 100 $REF $QRY && touch polish_align.success && rm -f filter.success || exit
fi

#delta-filter
if [ ! -e polish_filter.success ];then
parallel_delta-filter.sh $DELTAFILE '-r -l 100  -o 20 ' 9 && mv $DELTAFILE.fdelta $DELTAFILE.1.delta && \
touch polish_filter.success && rm -f  polish_add_not_aligning.success || exit
fi

if [ ! -e polish_add_not_aligning.success ];then
#add the sequences that did not align and longer than 1000 bp
show-coords -lcH -I $SIMILARITY_RATE $DELTAFILE.1.delta| perl -ane '{$palign{$F[-1]}+=$F[-4];}END{foreach $k(keys %palign){print $k,"\n" if($palign{$k}>20)}}' > aligned_sequences.txt
ufasta sizes -H $QRY | awk '{if($2<1000) print $1}' > short_sequences.txt
ufasta extract -v -f <(cat aligned_sequences.txt short_sequences.txt) $QRY > $REFN.$QRYN.extra.fa && \
cat $REF $REFN.$QRYN.extra.fa > $REFN.$QRYN.all.fa && \
touch polish_add_not_aligning.success && rm -f polish_replace_consensus.success || exit
fi

if [ ! -e polish_replace_consensus.success ];then
if [ $MERGE -gt 0 ];then
show-coords -lcHr -I $SIMILARITY_RATE -L 100 $DELTAFILE.1.delta |  merge_matches_and_tile_coords_file.pl $MERGE | reconcile_consensus.pl $REFN.$QRYN.all.fa $QRY > $REFN.$QRYN.all.polished.fa && \
touch polish_replace_consensus.success && rm -f polish_self_map.success || exit
else
show-coords -lcHr -I $SIMILARITY_RATE -L 100 $DELTAFILE.1.delta |  reconcile_consensus.pl $REFN.$QRYN.all.fa $QRY > $REFN.$QRYN.all.polished.fa && \
touch polish_replace_consensus.success && rm -f polish_self_map.success || exit
fi
fi

#here we map the contigs against themselves to figure out which ones are redundant
if [ ! -e polish_self_map.success ];then
perl -ane 'BEGIN{$seq=""}{if($F[0] =~ /^\>/){if(length($seq)>=1000){print length($seq)," $rn $seq\n";}$seq="";$rn=$F[0];}else{$seq.=$F[0];}}END{print length($seq)," $rn $seq\n";}' $REFN.$QRYN.all.polished.fa | sort -S 50% -nrk1 | awk '{print $2"\n"$3}' >  scaffolds.ref.fa && \
ufasta extract -f <(ufasta sizes -H $REFN.$QRYN.all.polished.fa | awk '{if($2<100000) print $1}') $REFN.$QRYN.all.polished.fa > $REFN.$QRYN.all.polished.lt100000.fa && \
nucmer -t $NUM_THREADS --batch $(($(stat -c%s "scaffolds.ref.fa")/5)) -c 100 -b 100 -p sasm_to_sasm  scaffolds.ref.fa  $REFN.$QRYN.all.polished.lt100000.fa && \
touch polish_self_map.success && rm -f polish_filter_map.success || exit
fi

if [ ! -e polish_filter_map.success ];then
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' sasm_to_sasm.delta > sasm_to_sasm.noself.delta && \
delta-filter -i $SIMILARITY_RATE -q -o 20 sasm_to_sasm.noself.delta > sasm_to_sasm.noself.fdelta && \
show-coords -lcHr sasm_to_sasm.noself.fdelta | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl 1000 | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' > sduplicates.txt && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' sasm_to_sasm.delta| show-coords -lcH /dev/stdin | awk '{if($12>$13 && $10>98 && $16>90) print $NF}' >> sduplicates.txt && \
ufasta extract -v -f sduplicates.txt $REFN.$QRYN.all.polished.fa > $REFN.$QRYN.all.polished.deduplicated.fa && \
touch polish_filter_map.success || exit
fi

echo "Output sequences in $REFN.$QRYN.all.polished.deduplicated.fa"
ufasta n50 -A -S -N50 $REFN.$QRYN.all.polished.deduplicated.fa
