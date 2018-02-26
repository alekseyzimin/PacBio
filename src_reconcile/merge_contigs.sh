#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1

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
        -r|--reference)
            REF="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: merge_contigs.sh -r <sequences to be merged> -q <sequence to merge with> -t <number of threads> "
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
if [ ! -e merge_align.success ];then
nucmer -t $NUM_THREADS -p  $DELTAFILE -l 31 -c 200 $REF $QRY && touch merge_align.success && rm -f merge_filter.success || exit
fi

#delta-filter
if [ ! -e merge_filter.success ];then
parallel_delta-filter.sh $DELTAFILE '-i 97 -r -l 1000' 9 && mv $DELTAFILE.fdelta $DELTAFILE.r.delta && \
touch merge_filter.success && rm -f  merge_merge.success || exit
fi

if [ ! -e merge_merge.success ];then
show-coords -lcHq -I 99 -L 1000 $DELTAFILE.r.delta | extract_merges.pl $QRY > merges.txt && merge_contigs.pl < merges.txt| create_merged_sequences.pl $REF merges.txt > $REFN.$QRYN.merged.fa && touch merge_merge.success || exit
fi

echo "Output sequences in $REFN.$QRYN.merged.fa"
ufasta n50 -A -S -N50 $REFN.$QRYN.merged.fa
