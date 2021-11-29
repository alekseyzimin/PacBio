#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
MIN_MATCH=5000
OVERHANG=1000
MIN_SCORE=60
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

function filter_convert_paf () {
#extract alignments of all long reads that satisfy the overhng, score and match length requirements and align to two or more contigs and convert to coords format
  awk '{ max_overhang=int("'$OVERHANG'"); if($4-$3>int("'$MIN_MATCH'") && ($8 < max_overhang || $7-$9 < max_overhang) && $12>=int("'$MIN_SCORE'")) print $0}' $1 |\
  perl -ane 'BEGIN{my %to_output=();}{
    push(@lines,join(" ",@F));if(not(defined($ctg{$F[0]}))){$ctg{$F[0]}=$F[5];}else{$to_output{$F[0]}=1 if(not($ctg{$F[0]} eq $F[5]));}
    }
    END{
    foreach $l(@lines){my @f=split(/\s+/,$l);print "$l\n" if(defined($to_output{$f[0]}));}
    }' | \
  sort -k1,1 -k3,3n -S 10% | \
  awk '{if($5=="+"){print $8+1" "$9" | "$3+1" "$4" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}else{print $8+1" "$9" | "$4" "$3+1" | "$11" "$4-$3" | 100 | "$7" "$2" | "int($11/$7*10000)/100" "int(($4-$3)/$2*10000)/100" | "$6" "$1}}' > $2.tmp && \
  mv $2.tmp $2
}

function usage () {
  echo "Usage:"
  echo "samba.sh -r <contigs or scaffolds in fasta format> -q <long reads or another assembly used to scaffold in fasta format> -t <number of threads> -m <minimum matching length, default:5000> -o <maximum overhang, default:1000> -a <optional: allowed merges file in the format per line: contig1 contig2, only pairs of contigs listed will be considered for merging, useful for intrascaffold gap filling>"
}

if [ $# -lt 1 ];then
  usage
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
            usage
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
$MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REF $QRY 1> $REFN.$QRYN.paf.tmp 2>minimap.err && mv $REFN.$QRYN.paf.tmp $REFN.$QRYN.paf && touch scaffold_align.success && rm -f scaffold_split.success || error_exit "minimap2 failed"
fi

if [ ! -e scaffold_split.success ];then
log "Filtering alignments and looking for misassemblies"
#first we figure out which reads we are keepping for subsequent steps.  We are keeping alignments of all reads that satisfy minimum alignment length criteria and map to 2+ contigs
awk '{if($4-$3>int("'$MIN_MATCH'") && $12>=int("'$MIN_SCORE'")) print $0}' $REFN.$QRYN.paf | \
perl -ane '{push(@lines,join(" ",@F));$ctg{$F[0]}.="$F[5] ";}
END{
my %to_output=();
foreach $r (keys %ctg){my @f=split(/\s+/,$ctg{$r}); 
my %temp=(); 
foreach $c(@f){$temp{$c}=1}; 
my $size = keys %temp; if($size>1){$to_output{$r}=1;}
}
foreach $l(@lines){my @f=split(/\s+/,$l);print "$l\n" if(defined($to_output{$f[0]}));}}' | \
sort -k1,1 -k3,3n -S 10% > $REFN.$QRYN.filtered.paf.tmp && mv $REFN.$QRYN.filtered.paf.tmp $REFN.$QRYN.filtered.paf && \
awk 'BEGIN{r="";c="";oh=int("'$MIN_MATCH'");}{
if($1==r && $6!=c){
split(prevline,a," ");
if(a[5]=="+"){if(a[7]-a[9] >= oh && a[2]-a[4] >= oh){print a[6]" "a[9];}}else{ if(a[8] >= oh && a[2]-a[4] >= oh){print a[6]" "a[8]}} 
  if($5=="+"){       if($8 >= oh &&    $3 >= oh){print $6" "$8;}}    else{if($7-$9 >= oh &&    $3 >= oh){print $6" "$9;}}
}prevline=$0;c=$6;r=$1;}' $REFN.$QRYN.filtered.paf | \
sort -S10% -k1,1 -k2,2n | \
uniq -c |awk '{if($1<3) print $2" "$3}' | \
perl -ane '{push(@ctg,$F[0]);push(@coord,$F[1]);}END{$rad=30;$tol=500;for($i=0;$i<=$#ctg;$i++){my $score=0;for($j=$i-$rad;$j<$i+$rad;$j++){next if($j<0 || $j>$#ctg);if(abs($coord[$j]-$coord[$i])<$tol && $ctg[$j] eq $ctg[$i]){$score+=exp(-abs($coord[$j]-$coord[$i])/$tol)}} print "$ctg[$i] $coord[$i] $score\n" if($score>3);}}' | \
sort -nrk3,3 -S10% | uniq | perl -ane '{if(not(defined($h{$F[0]}))){$h{$F[0]}=$F[1];}else{@f=split(/\s+/,$h{$F[0]});my $flag=0;foreach $v(@f){$flag=1 if(abs($v-$F[1])<int("'$MIN_MATCH'"));}if($flag==0){$h{$F[0]}.=" $F[1]"}}}END{foreach $k(keys %h){@f=split(/\s+/,$h{$k});foreach $v(@f){print "break $k $v\n"}}}' |
sort -S10% -k2,2 -k3,3n > $REFN.$QRYN.splits.txt && \
echo -n "Found misassemblies: " && wc -l $REFN.$QRYN.splits.txt && \
$MYPATH/break_contigs.pl  $REFN.$QRYN.splits.txt < $REF > $REFN.split.fa.tmp && mv $REFN.split.fa.tmp $REFN.split.fa && touch scaffold_split.success && \
rm -f scaffold_split_align.success || error_exit "splitting at misassemblies failed" 
fi

#minimap
if [ ! -e scaffold_split_align.success ];then
log "Aligning the reads to the split contigs"
$MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REFN.split.fa <(ufasta extract -f <(awk '{print $1}' $REFN.$QRYN.filtered.paf) $QRY) 1> $REFN.$QRYN.split.paf.tmp 2>minimap.err && mv $REFN.$QRYN.split.paf.tmp $REFN.$QRYN.split.paf && touch scaffold_split_align.success && rm -f scaffold_filter.success || error_exit "minimap2 failed"
fi

if [ ! -e scaffold_filter.success ];then
filter_convert_paf $REFN.$QRYN.split.paf $REFN.$QRYN.coords || error_exit "filtering alignments failed"
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
$MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -x map-ont $REFN.split.fa $REFN.$QRYN.patches.fa 2>minimap.err  > $REFN.$QRYN.patches.paf.tmp && mv  $REFN.$QRYN.patches.paf.tmp  $REFN.$QRYN.patches.paf && touch scaffold_align_patches.success && rm -f scaffold_scaffold.success || error_exit "minimap2 patches failed"
fi

if [ ! -e scaffold_scaffold.success ];then
log "Creating scaffold graph and building scaffolds"
rm -f do_consensus.sh && \
filter_convert_paf $REFN.$QRYN.patches.paf $REFN.$QRYN.patches.coords && \
cat $REFN.$QRYN.patches.coords | $MYPATH/extract_merges.pl $REFN.$QRYN.patches.fa $ALLOWED > $REFN.$QRYN.patches.links.txt.tmp && \
mv $REFN.$QRYN.patches.links.txt.tmp $REFN.$QRYN.patches.links.txt && \
$MYPATH/find_repeats.pl $REFN.$QRYN.patches.coords $REFN.$QRYN.patches.links.txt >$REFN.repeats.txt.tmp && \
mv $REFN.repeats.txt.tmp $REFN.repeats.txt && \
perl -ane '$h{$F[0]}=1;END{open(FILE,"'$REFN.$QRYN.patches.coords'");while($line=<FILE>){@f=split(/\s+/,$line);print $line unless(defined($h{$f[-2]}));}}' $REFN.repeats.txt | \
$MYPATH/extract_merges.pl $REFN.$QRYN.patches.fa $ALLOWED > $REFN.$QRYN.patches.uniq.links.txt.tmp && \
mv $REFN.$QRYN.patches.uniq.links.txt.tmp $REFN.$QRYN.patches.uniq.links.txt && \
$MYPATH/merge_contigs.pl $REFN.split.fa < $REFN.$QRYN.patches.uniq.links.txt 2>$REFN.$QRYN.bubbles.txt | \
$MYPATH/insert_repeats.pl $REFN.repeats.txt |\
$MYPATH/create_merged_sequences.pl $REFN.split.fa  <(cat $REFN.$QRYN.patches.uniq.links.txt $REFN.$QRYN.patches.links.txt |sort -S 10% |uniq) | \
$MYPATH/ufasta extract -v -f $REFN.$QRYN.bubbles.txt > $REFN.$QRYN.scaffolds.fa.tmp && \
mv $REFN.$QRYN.scaffolds.fa.tmp $REFN.scaffolds.all.fa  && touch scaffold_scaffold.success && rm -f scaffold_rejoin.success || error_exit "walking the scaffold graph failed"
fi

if [ ! -e scaffold_rejoin.success ];then
ufasta sizes -H $REFN.scaffolds.all.fa | $MYPATH/make_rejoin_links.pl > $REFN.rejoin.links.txt.tmp && mv $REFN.rejoin.links.txt.tmp $REFN.rejoin.links.txt && \
awk '{print $1" "$2" "$3" "$4" "$5" "$6}' $REFN.rejoin.links.txt | \
$MYPATH/create_merged_sequences.pl $REFN.scaffolds.all.fa $REFN.rejoin.links.txt > $REFN.scaffolds.all.rejoin.fa.tmp && mv $REFN.scaffolds.all.rejoin.fa.tmp $REFN.scaffolds.all.rejoin.fa && touch scaffold_rejoin.success && rm -f scaffold_deduplicate.success || error_exit "rejoining spurios breaks failed"
fi

if [ ! -e scaffold_deduplicate.success ];then
log "Deduplicating contigs"
awk 'BEGIN{n=0}{if(length($NF)>100){print ">"n"\n"$NF;n++}}' $REFN.$QRYN.patches.uniq.links.txt > $REFN.$QRYN.patches.uniq.links.fa.tmp && mv $REFN.$QRYN.patches.uniq.links.fa.tmp $REFN.$QRYN.patches.uniq.links.fa && \
MAX_PATCH=`ufasta sizes -H  $REFN.$QRYN.patches.uniq.links.fa | sort -nrk2,2 -S 10% | head -n 1 | awk '{print $2}'`
$MYPATH/ufasta extract -f <($MYPATH/ufasta sizes -H $REFN.scaffolds.all.rejoin.fa | awk '{if($2<int("'$MAX_PATCH'")) print $1}') $REFN.scaffolds.all.rejoin.fa > $REFN.short_contigs.fa.tmp && mv $REFN.short_contigs.fa.tmp $REFN.short_contigs.fa && \
ufasta extract -v -f <($MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS $REFN.short_contigs.fa $REFN.$QRYN.patches.uniq.links.fa 2>/dev/null | awk '{if(($9-$8)/$7>.90) print $6}') $REFN.scaffolds.all.rejoin.fa > $REFN.scaffolds.fa.tmp && mv $REFN.scaffolds.fa.tmp $REFN.scaffolds.fa && \
rm $REFN.short_contigs.fa $REFN.$QRYN.patches.uniq.links.fa && touch scaffold_deduplicate.success || error_exit "deduplicate failed"
fi

log "Output scaffold sequences in $REFN.scaffolds.fa"
