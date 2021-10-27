#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
MIN_MATCH=5000
OVERHANG=1000
ALLOWED=""
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

if [ $# -lt 1 ];then
  error_exit "Usage: samba.sh -r <contigs or scaffolds in fasta format> -q <long reads or another assembly used to scaffold in fasta format> -t <number of threads> -m <minimum matching length, default:5000> -o <maximum overhang, default:1000> -a <optional: allowed merges file in the format per line: contig1 contig2, only pairs of contigs listed will be considered for merging, useful for intrascaffold gap filling>"
fi

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
        -a|--allowed)
            ALLOWED="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: samba.sh -r <contigs or scaffolds in fasta format> -q <long reads or another assembly used to scaffold in fasta format> -t <number of threads> -m <minimum matching length, default:5000> -o <maximum overhang, default:1000> -a <optional: allowed merges file in the format per line: contig1 contig2, only pairs of contigs listed will be considered for merging, useful for intrascaffold gap filling>"
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
awk '{max_overhang=0.1*$7;if(max_overhang>int("'$OVERHANG'")) max_overhang=int("'$OVERHANG'"); if($4-$3>int("'$MIN_MATCH'") && ($8<max_overhang || $7-$9<max_overhang) && $12>=60) print $0}' $REFN.$QRYN.paf | \
sort -k1,1 -k3,3n -S 10% | \
awk 'BEGIN{r="";c=""}{if($1==r && $6!=c){print prevline"\n"$0;}prevline=$0;c=$6;r=$1;}' | \
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
cat  $REFN.$QRYN.coords |$MYPATH/extract_merges.pl $REFN.$QRYN.reads.fa $ALLOWED > $REFN.$QRYN.links.txt.tmp && mv $REFN.$QRYN.links.txt.tmp $REFN.$QRYN.links.txt && \
$MYPATH/find_repeats.pl $REFN.$QRYN.coords $REFN.$QRYN.links.txt >$REFN.repeats.txt.tmp && mv $REFN.repeats.txt.tmp $REFN.repeats.txt && \
rm -f patches.polished.fa && \
echo "#!/bin/bash" > do_consensus.sh && \
echo "rm -rf polish.?.tmp" >> do_consensus.sh
for jindex in $(seq 0 9);do
echo "if [ -s patches.ref.$jindex.fa ] && [ -s patches.reads.$jindex.fa ];then $MYPATH/../Flye/bin/flye --polish-target patches.ref.$jindex.fa --iterations 1 --nano-raw patches.reads.$jindex.fa --threads 8 --out-dir polish.$jindex.tmp 1>polish.err 2>&1 && rm -f patches.ref.$jindex.fa patches.reads.$jindex.fa;fi &"  >> do_consensus.sh 
done
echo "wait;" >> do_consensus.sh && \
echo "cat polish.?.tmp/polished_1.fasta | $MYPATH/ufasta one >> patches.polished.fa" >> do_consensus.sh && \
chmod 0755 do_consensus.sh && \
perl -ane '$h{$F[0]}=1;END{open(FILE,"'$REFN.$QRYN.coords'");while($line=<FILE>){@f=split(/\s+/,$line);print $line unless(defined($h{$f[-2]}));}}' $REFN.repeats.txt | \
$MYPATH/extract_merges.pl $REFN.$QRYN.reads.fa  $ALLOWED >/dev/null && \
rm -f do_consensus.sh && \
touch patches.polished.fa && \
cat patches.polished.fa patches.raw.fa > $REFN.$QRYN.patches.fa.tmp && mv $REFN.$QRYN.patches.fa.tmp $REFN.$QRYN.patches.fa && \
rm -rf patches.ref.fa patches.reads.fa patches.raw.fa polish.?.tmp && \
touch scaffold_links.success && rm -f scaffold_align_patches.success || error_exit "links consensus failed"
fi

if [ ! -e scaffold_align_patches.success ];then
log "Aligning the scaffolding sequences to the contigs"
$MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REF $REFN.$QRYN.patches.fa 2>minimap.err  > $REFN.$QRYN.patches.paf.tmp && mv  $REFN.$QRYN.patches.paf.tmp  $REFN.$QRYN.patches.paf && touch scaffold_align_patches.success && rm -f scaffold_scaffold.success || error_exit "minimap2 patches failed"
fi

if [ ! -e scaffold_scaffold.success ];then
log "Creating scaffold graph and building scaffolds"
rm -f do_consensus.sh && \
awk '{if($4-$3>int("'$MIN_MATCH'") && ($8<int("'$OVERHANG'") || $7-$9<int("'$OVERHANG'")) && $12>=60) print $0}' $REFN.$QRYN.patches.paf |\
sort -k1,1 -k3,3n -S 10% | \
awk 'BEGIN{r="";c=""}{if($1==r && $6!=c){print prevline"\n"$0;}prevline=$0;c=$6;r=$1;}' | \
awk '{if($5=="+"){print $8+1" "$9" | "$3+1" "$4" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}else{print $8+1" "$9" | "$4" "$3+1" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}}' |\
$MYPATH/extract_merges.pl $REFN.$QRYN.patches.fa $ALLOWED > $REFN.$QRYN.patches.links.txt.tmp && mv $REFN.$QRYN.patches.links.txt.tmp $REFN.$QRYN.patches.links.txt && \
$MYPATH/find_repeats.pl $REFN.$QRYN.patches.coords $REFN.$QRYN.patches.links.txt >$REFN.repeats.txt.tmp && mv $REFN.repeats.txt.tmp $REFN.repeats.txt && \
perl -ane '$h{$F[0]}=1;END{open(FILE,"'$REFN.$QRYN.patches.coords'");while($line=<FILE>){@f=split(/\s+/,$line);print $line unless(defined($h{$f[-2]}));}}' $REFN.repeats.txt | \
$MYPATH/extract_merges.pl $REFN.$QRYN.patches.fa $ALLOWED > $REFN.$QRYN.patches.uniq.links.txt.tmp && mv $REFN.$QRYN.patches.uniq.links.txt.tmp $REFN.$QRYN.patches.uniq.links.txt && \
$MYPATH/merge_contigs.pl $REF < $REFN.$QRYN.patches.uniq.links.txt 2>$REFN.$QRYN.bubbles.txt | \
$MYPATH/insert_repeats.pl $REFN.repeats.txt |\
$MYPATH/create_merged_sequences.pl $REF  <(cat $REFN.$QRYN.patches.uniq.links.txt $REFN.$QRYN.patches.links.txt |sort -S 10% |uniq) | \
$MYPATH/ufasta extract -v -f $REFN.$QRYN.bubbles.txt > $REFN.$QRYN.scaffolds.fa.tmp && \
mv $REFN.$QRYN.scaffolds.fa.tmp $REFN.scaffolds.all.fa  && touch scaffold_scaffold.success && rm -f scaffold_deduplicate.success || error_exit "walking the scaffold graph failed"
fi

if [ ! -e scaffold_deduplicate.success ];then
log "Deduplicating contigs"
awk 'BEGIN{n=0}{if(length($NF)>100){print ">"n"\n"$NF;n++}}' $REFN.$QRYN.patches.uniq.links.txt > $REFN.$QRYN.patches.uniq.links.fa.tmp && mv $REFN.$QRYN.patches.uniq.links.fa.tmp $REFN.$QRYN.patches.uniq.links.fa && \
$MYPATH/ufasta extract -f <($MYPATH/ufasta sizes -H $REFN.scaffolds.all.fa | awk '{if($2<int("'$MIN_MATCH'")) print $1}') $REFN.scaffolds.all.fa > $REFN.short_contigs.fa.tmp && mv $REFN.short_contigs.fa.tmp $REFN.short_contigs.fa && \
ufasta extract -v -f <($MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS $REFN.short_contigs.fa $REFN.$QRYN.patches.uniq.links.fa 2>/dev/null | awk '{if(($9-$8)/$7>.95) print $6}') $REFN.scaffolds.all.fa > $REFN.scaffolds.fa.tmp && mv $REFN.scaffolds.fa.tmp $REFN.scaffolds.fa && \
rm $REFN.short_contigs.fa $REFN.$QRYN.patches.uniq.links.fa && touch scaffold_deduplicate.success || error_exit "deduplicate failed"
fi

log "Output scaffold sequences in $REFN.scaffolds.fa"