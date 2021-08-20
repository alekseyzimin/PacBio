#!/bin/bash
#this script closes gaps in chromosome scaffolds using another assembly or the reference
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
MIN_MATCH=2500
OVERHANG=1000
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
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
        -i|--identity)
            IDENTITY="$2"
            shift
            ;;
        -q|--query)
            QRY="$2"
            shift
            ;;
        -m|--min-match)
            MIN_MATCH="$2"
            shift
            ;;
        -r|--reference)
            REF="$2"
            shift
            ;;
        -o|--overhang)
            OVERHANG="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: close_scaffold_gaps.sh <options>"
            echo "-r <scaffolds to gapclose> MANDATORY"
            echo "-q <gapclosing sequences, can be long reads> MANDATORY"
            echo "-t <number of threads> default:1"
            echo "-i <identity%> default:98"
            echo "-m <minimum match length on the two sides of the gap> default:2500"
            echo "-o <max overhang> default:1000"
            echo "-v verbose"
            echo "-h|--help|-u|--usage this message"
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
error_exit "merging sequence $QRY does not exist or size zero"
fi

REFN=`basename $REF`
QRYN=`basename $QRY`

#split
if [ ! -e "scaffold_merge.split.success" ];then
log "Splitting scaffolds into contigs"
$MYPATH/splitScaffoldsAtNs.pl  < $REF | ufasta one > $REFN.split && \
grep '^>' --text $REFN.split | perl -ane '{($rn,$coord)=split(/\./,substr($F[0],1));$h{$rn}.=substr($F[0],1)." ";}END{foreach $r(keys %h){@f=split(/\s+/,$h{$r}); for ($i=0;$i<$#f;$i++){print $f[$i]," ",$f[$i+1],"\n"}}}' > $REFN.valid_join_pairs.txt && \
touch scaffold_merge.split.success && rm -rf scaffold_merge.scaffold.success
fi

if [ ! -e "scaffold_merge.scaffold.success" ];then
log "Closing gaps"
$MYPATH/masurca_scaffold.sh -r $REFN.split -q $QRY -t $NUM_THREADS -o $OVERHANG -m $MIN_MATCH -a $REFN.valid_join_pairs.txt && \
$MYPATH/recover_scaffolds.pl < $REFN.split.scaffolds.fa |ufasta format > $REFN.split.joined.tmp && \
mv $REFN.split.joined.tmp $REFN.split.joined.fa && \
touch scaffold_merge.scaffold.success
fi

log "Output sequences in $REFN.split.joined.fa"
