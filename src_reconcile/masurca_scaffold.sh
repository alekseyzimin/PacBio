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
            echo "Usage: masurca_scaffold.sh -r <sequences to be merged> -q <sequence to merge with> -t <number of threads> -m <minimum matching length, default:5000> -o <maximum overhang, default:1000>"
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
if [ ! -e scaffold_align.success ];then
log "Aligning the reads to the contigs"
$MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REF $QRY 1> $REFN.$QRYN.paf.tmp 2>minimap.err && mv $REFN.$QRYN.paf.tmp $REFN.$QRYN.paf && touch scaffold_align.success && rm -f scaffold_filter.success || error_exit "minimap2 failed"
fi

if [ ! -e scaffold_filter.success ];then
log "Filtering alignments"
awk '{if($4-$3>int("'$MIN_MATCH'") && ($8<int("'$OVERHANG'") || $7-$9<int("'$OVERHANG'")) && $12>=60) print $0}' $REFN.$QRYN.paf | \
sort -k1,1 -k3,3n -S 10% | \
awk 'BEGIN{r="";c=""}{if($1!=r){print $0" "$1;r=$1;c=$6}else if($6!=c){print $0" "$1;c=$6}}' | \
uniq -D -f 18 | \
awk '{if($5=="+"){print $8+1" "$9" | "$3+1" "$4" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}else{print $8+1" "$9" | "$4" "$3+1" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}}' > $REFN.$QRYN.coords.tmp && mv $REFN.$QRYN.coords.tmp $REFN.$QRYN.coords || error_exit "filtering alignments failed" 
touch scaffold_filter.success && rm -f  scaffold_reads.success 
fi

if [ ! -e scaffold_reads.success ];then
log "Extracting reads for the patches"
$MYPATH/ufasta extract -f <(awk '{print $NF}' $REFN.$QRYN.coords) $QRY > $REFN.$QRYN.reads.fa.tmp && mv $REFN.$QRYN.reads.fa.tmp $REFN.$QRYN.reads.fa && \
touch scaffold_reads.success && rm -f  scaffold_links.success || error_exit "failed in extracting the reads for scaffolding"
fi

if [ ! -e scaffold_links.success ];then
#this is a bit tricky
#we first find all links to identify repeats by both coverage and linking criteria
#then we exclude the repeats and re-compute the links to create input for consensus patches
#the reasoning is that we do not want to create extra patches for repetitive junctions
#we now create do_consensus.sh and extract_merges.pl will see it and will use it do the consensus for the patches, do not forgrt to delete it
log "Creating scaffold links"
rm -f do_consensus.sh && \
cat  $REFN.$QRYN.coords |$MYPATH/extract_merges.pl $REFN.$QRYN.reads.fa  > $REFN.$QRYN.links.txt.tmp && mv $REFN.$QRYN.links.txt.tmp $REFN.$QRYN.links.txt && \
$MYPATH/find_repeats.pl $REFN.$QRYN.coords $REFN.$QRYN.links.txt >$REFN.repeats.txt.tmp && mv $REFN.repeats.txt.tmp $REFN.repeats.txt && \
rm -f patches.polished.fa && \
echo "#!/bin/bash" > do_consensus.sh && \
echo "rm -rf polish.tmp" >> do_consensus.sh && \
echo "$MYPATH/../Flye/bin/flye --polish-target patches.ref.fa --iterations 1 --nano-raw patches.reads.fa --threads $NUM_THREADS --out-dir polish.tmp 1>polish.err 2>&1 && $MYPATH/ufasta one polish.tmp/polished_1.fasta >> patches.polished.fa"  >> do_consensus.sh && \
chmod 0755 do_consensus.sh && \
perl -ane '$h{$F[0]}=1;END{open(FILE,"'$REFN.$QRYN.coords'");while($line=<FILE>){@f=split(/\s+/,$line);print $line unless(defined($h{$f[-2]}));}}' $REFN.repeats.txt | \
$MYPATH/extract_merges.pl $REFN.$QRYN.reads.fa  >/dev/null && \
rm -f do_consensus.sh && \
cat patches.polished.fa patches.raw.fa > $REFN.$QRYN.patches.fa.tmp && mv $REFN.$QRYN.patches.fa.tmp $REFN.$QRYN.patches.fa && \
rm -rf patches.ref.fa patches.reads.fa patches.raw.fa polish.tmp && \
touch scaffold_links.success && rm -f scaffold_align_patches.success || error_exit "links consensus failed"
fi

if [ ! -e scaffold_align_patches.success ];then
log "Aligning the scaffolding sequences to the contigs"
$MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REF $REFN.$QRYN.patches.fa 2>minimap.err | \
awk '{if($4-$3>int("'$MIN_MATCH'") && ($8<int("'$OVERHANG'") || $7-$9<int("'$OVERHANG'")) && $12>=60) print $0}' |\
sort -k1,1 -k3,3n -S 10% | \
awk 'BEGIN{r="";c=""}{if($1!=r){print $0" "$1;r=$1;c=$6}else if($6!=c){print $0" "$1;c=$6}}' | \
uniq -D -f 18 | \
awk '{if($5=="+"){print $8+1" "$9" | "$3+1" "$4" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}else{print $8+1" "$9" | "$4" "$3+1" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}}' > $REFN.$QRYN.patches.coords.tmp && mv $REFN.$QRYN.patches.coords.tmp $REFN.$QRYN.patches.coords && \
touch scaffold_align_patches.success || error_exit "minimap2 patches failed"
fi

if [ -e scaffold_align_patches.success ];then
log "Creating scaffold graph and building scaffolds"
rm -f do_consensus.sh && \
cat $REFN.$QRYN.patches.coords | $MYPATH/extract_merges.pl $REFN.$QRYN.patches.fa > $REFN.$QRYN.patches.links.txt.tmp && mv $REFN.$QRYN.patches.links.txt.tmp $REFN.$QRYN.patches.links.txt && \
$MYPATH/find_repeats.pl $REFN.$QRYN.coords $REFN.$QRYN.patches.links.txt >$REFN.repeats.txt.tmp && mv $REFN.repeats.txt.tmp $REFN.repeats.txt && \
perl -ane '$h{$F[0]}=1;END{open(FILE,"'$REFN.$QRYN.patches.coords'");while($line=<FILE>){@f=split(/\s+/,$line);print $line unless(defined($h{$f[-2]}));}}' $REFN.repeats.txt | \
$MYPATH/extract_merges.pl $REFN.$QRYN.patches.fa > $REFN.$QRYN.patches.uniq.links.txt.tmp && mv $REFN.$QRYN.patches.uniq.links.txt.tmp $REFN.$QRYN.patches.uniq.links.txt && \
$MYPATH/merge_contigs.pl $REF < $REFN.$QRYN.patches.uniq.links.txt 2>$REFN.$QRYN.bubbles.txt | \
$MYPATH/insert_repeats.pl $REFN.repeats.txt |\
$MYPATH/create_merged_sequences.pl $REF  <(cat $REFN.$QRYN.patches.uniq.links.txt $REFN.$QRYN.patches.links.txt |sort -S 10% |uniq) | \
$MYPATH/ufasta extract -v -f $REFN.$QRYN.bubbles.txt > $REFN.$QRYN.scaffolds.fa.tmp && \
mv $REFN.$QRYN.scaffolds.fa.tmp $REFN.scaffolds.fa || error_exit "walking the scaffold graph failed"
fi

log "Output scaffold sequences in $REFN.scaffolds.fa"
