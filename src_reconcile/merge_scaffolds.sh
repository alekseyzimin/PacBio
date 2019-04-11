#!/bin/bash
#this script closes gaps in chromosome scaffolds using another assembly or the reference
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
MIN_MATCH=1000
IDENTITY=98
OVERHANG=5000
MINGAP=1000
MAXGAP=50000

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
        -G|--maxgap)
            MAXGAP="$2"
            shift
            ;;
        -g|--mingap)
            MINGAP=$((-1*$2))
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: merge_scaffolds.sh -r <sequences to be merged> -q <sequence to merge with> -t <number of threads> "
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

#split
if [ ! -e "scaffold_merge.split.success" ];then
$MYPATH/splitScaffoldsAtNs.pl  < $REF > $REFN.split && \
grep '^>' --text $REFN.split | perl -ane '{($rn,$coord)=split(/\./,substr($F[0],1));$h{$rn}.=substr($F[0],1)." ";}END{foreach $r(keys %h){@f=split(/\s+/,$h{$r}); for ($i=0;$i<$#f;$i++){print $f[$i]," ",$f[$i+1],"\n"}}}' > valid_join_pairs.txt && \
touch scaffold_merge.split.success && rm -rf scaffold_merge.align.success
fi

#nucmer
if [ ! -e scaffold_merge.align.success ];then
nucmer -t $NUM_THREADS -p  $DELTAFILE -l 31 -c 200 $REFN.split $QRY && \
touch scaffold_merge.align.success && \
rm -f scaffold_merge.filter.success || exit
fi

#delta-filter
if [ ! -e scaffold_merge.filter.success ];then
parallel_delta-filter.sh $DELTAFILE "-1 -l 200" 9 && mv $DELTAFILE.fdelta $DELTAFILE.r.delta && \
touch scaffold_merge.filter.success && rm -f  scaffold_merge.merge.success || exit
fi

#perform merge
if [ ! -e scaffold_merge.merge.success ];then
$MYPATH/show-coords -lcHq -I $IDENTITY $DELTAFILE.r.delta | $MYPATH/extract_merges_mega-reads.pl $QRY  valid_join_pairs.txt $OVERHANG $MINGAP $MAXGAP > merges.txt && \
perl -ane '{if($F[2] eq "F"){$merge="$F[0] $F[3]";}else{$merge="$F[3] $F[0]";} if(not(defined($h{$merge}))|| $h{$merge} > $F[1]+$F[4]){$hl{$merge}=join(" ",@F);$h{$merge}=$F[1]+$F[4];}}END{foreach $k(keys %hl){print $hl{$k},"\n"}}' merges.txt > merges.best.txt && \
cat \
<($MYPATH/ufasta extract -v -f <(awk '{print $1"\n"$2;}' valid_join_pairs.txt) $REFN.split) \
<($MYPATH/merge_mega-reads.pl < merges.best.txt | $MYPATH/create_merged_mega-reads.pl <($MYPATH/ufasta extract -f <(awk '{print $1"\n"$2;}' valid_join_pairs.txt) $REFN.split) merges.best.txt) | \
$MYPATH/recover_scaffolds.pl > $REFN.split.joined.tmp && \
mv $REFN.split.joined.tmp $REFN.split.joined.fa && \
touch scaffold_merge.merge.success
fi

echo "Output sequences in $REFN.split.joined.fa"
