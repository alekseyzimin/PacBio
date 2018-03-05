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
        -m|--merge)
            MERGE_SLACK="$2"
            shift
            ;;
        -p|--path)
            SMRT_PATH="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: polish_with_illumina_assembly_pbdagcon.sh -r <sequence to be polished> -q <sequence to polish with> -t <number of threads> -m <merge polishing sequence alignments slack (advanced)> "
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

if [ ! -e $SMRT_PATH/blasr ];then
error_exit "blasr not found at $SMRT_PATH"
fi


REFN=`basename $REF`
QRYN=`basename $QRY`
DELTAFILE=$REFN.$QRYN



#nucmer
if [ ! -e polish_align.success ];then
blasr  $QRY $REF -minMatch 15 -bestn 1 -m 5 -nproc $NUM_THREADS > $REFN.$QRYN.m5 2>/dev/null && touch polish_align.success || error_exit "blasr failed on the input data"
fi

if [ ! -e polish_replace_consensus.success ];then
$SMRT_PATH/pbdagcon -t 0 -m 100 -j $NUM_THREADS -c 1 $REFN.$QRYN.m5 > $REFN.$QRYN.cns && touch polish_replace_consensus.success || error_exit "pbdagcon failed on the input data"
fi

echo "Output sequences in $REFN.$QRYN.all.polished.deduplicated.fa"
ufasta n50 -A -S -N50 $REFN.$QRYN.all.polished.deduplicated.fa
