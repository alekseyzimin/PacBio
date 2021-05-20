#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
MIN_MATCH=5000
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
        -o|--overhang)
            OVERHANG="$2"
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
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: merge_contigs.sh -r <sequences to be merged> -q <sequence to merge with> -t <number of threads> -m <minimum matching length, default:1000> -i <matching identity, default:98>"
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
DELTAFILE=$REFN.$QRYN

#minimap
log "Aligning the scaffolding sequences to the contigs"
if [ ! -e scaffold_align.success ];then
$MYPATH/../FLye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REF $QRY 1> $REFN.$QRYN.paf.tmp 2>minimap.err && mv $REFN.$QRYN.paf.tmp $REFN.$QRYN.paf && touch scaffold_align.success && rm -f scaffold_filter.success || error_exit "minimap2 failed"
fi

#filter
log "Filtering alignments"
if [ ! -e scaffold_filter.success ];then
awk '{if($4-$3>int("'$MIN_MATCH'") && ($8<int("'$OVERHANG'") || $7-$9<int("'$OVERHANG'")) && $12>=60) print $0}' $REFN.$QRYN.paf | \
sort -k1,1 -k3,3n -S 10% | \
awk 'BEGIN{r="";c=""}{if($1!=r){print $0" "$1;r=$1;c=$6}else if($6!=c){print $0" "$1;c=$6}}' | \
uniq -D -f 18 | \
awk '{if($5=="+"){print $8+1" "$9+1" | "$3+1" "$4+1" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}else{print $8+1" "$9+1" | "$4+1" "$3+1" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}}' > $REFN.$QRYN.coords.tmp && mv $REFN.$QRYN.coords.tmp $REFN.$QRYN.coords || error_exit "filtering alignments failed" 
touch scaffold_filter.success && rm -f  scaffold_links.success 
fi

log "Creating scaffold links"
if [ ! -e scaffold_links.success ];then
cat  $REFN.$QRYN.coords |$MYPATH/extract_merges.pl <(ufasta extract -f <(awk '{print $NF}' $REFN.$QRYN.coords) $QRY) | \
awk '{print $7" "$8" "$2" "$5" "$2+$5" "$1" "$3" "$4" "$6}' | sort -k6,9 -k5,5n -S 10% | uniq -f 5 | awk '{print $6" "$3" "$7" "$8" "$4" "$9" "$1" "$2}'  > $REFN.$QRYN.links.txt.tmp && mv  $REFN.$QRYN.links.txt.tmp  $REFN.$QRYN.links.txt || error_exit "creating scaffold links failed"
touch scaffold_links.success && rm -f  scaffold_graph.success
fi

log "Creating scaffold graph and building scaffolds"
if [ ! -e scaffold_graph.success ];then
cat $REFN.$QRYN.links.txt | $MYPATH/merge_contigs.pl | $MYPATH/create_merged_sequences.pl $REF  $REFN.$QRYN.links.txt > $REFN.$QRYN.scaffolds.fa.tmp && mv $REFN.$QRYN.scaffolds.fa.tmp $REFN.$QRYN.scaffolds.fa || error_exit "walking the scaffold graph failed"
fi

log "Output sequences in $REFN.$QRYN.scaffolds.fa"
ufasta n50 -a $REFN.$QRYN.scaffolds.fa
