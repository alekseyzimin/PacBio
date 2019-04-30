#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
FLYE_PATH="$MYPATH/../Flye/bin";
ESTIMATED_GENOME_SIZE=0
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
MER=17
B=20
d=0.05
KMER=41
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
GC='\e[0;32m'
RC='\e[0;31m'
NC='\e[0m'
fi

log () {
  dddd=$(date)
  echo -e "${GC}[$dddd]${NC} $@"
}



function error_exit {
  echo "$1" >&2
    exit "${2:-1}"
}


#parsing arguments
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -m|--masurca_run_path)
    MASURCA_ASSEMBLY_WORK1_PATH="$2"
    shift
    ;;
    -t|--threads)
    NUM_THREADS="$2"
    shift
    ;;
    -M|--alignment_mer)
    MER="$2"
    shift
    ;;
    -B|--alignment_threshold)
    B="$2"
    shift
    ;;
    -D|--density)
    d="$2"
    shift
    ;;
    -r|--reference)
    REF="$2"
    shift
    ;;
    -a|--assembler_path)
    CA_PATH="$2"
    shift
    ;;
    -e|--estimated-genome-size)
    ESTIMATED_GENOME_SIZE="$2"
    shift
    ;;
    -o|--other_frg)
    OTHER_FRG="$2"
    shift
    ;;
    -v|--verbose)
    set -x
    ;;
    -h|--help|-u|--usage)
    echo "Usage: mega_reads_assemblei_ref.sh -m <path to MaSuRCA run work1 folder contents> -r <reference assembly fasta> -a <path to the location of runCA in wgs-8.2 instalation>"
    exit 0
    ;;
    *)
    echo "Unknown option $1"
    exit 1        # unknown option
    ;;
esac
shift
done

###############checking arguments#########################
if [ ! -e $CA_PATH/runCA ];then
echo "runCA not found at $CA_PATH!";
exit 1;
fi

export PATH=$CA_PATH:$MYPATH:$PATH:$FLYE_PATH

if [ ! -e $REF ];then
echo "Reference sequence file $REF not found!";
exit 1 ;
fi

COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}

log "Running synteny-guided assembly"
log "Using mer size $MER for mapping, B=$B, d=$d"
log "Using CA installation from $CA_PATH"
log "Using $NUM_THREADS threads"
log "Output prefix $COORDS"

rm -f .rerun

#first we re-create k-unitigs and super reads with smaller K
if [ ! -e superReadSequences.named.fasta ];then
log "Reducing super-read k-mer size"
awk 'BEGIN{n=0}{if($1~/^>/){}else{print ">sr"n"\n"$0;n+=2;}}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta > superReadSequences.fasta.in
create_k_unitigs_large_k -q 1 -c $(($KMER-1)) -t $NUM_THREADS -m $KMER -n $(($ESTIMATED_GENOME_SIZE*2)) -l $KMER -f `perl -e 'print 1/'$KMER'/1e5'` $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta.all  | grep --text -v '^>' | perl -ane '{$seq=$F[0]; $F[0]=~tr/ACTGactg/TGACtgac/;$revseq=reverse($F[0]); $h{($seq ge $revseq)?$seq:$revseq}=1;}END{$n=0;foreach $k(keys %h){print ">",$n++," length:",length($k),"\n$k\n"}}' > guillaumeKUnitigsAtLeast32bases_all.fasta.tmp && mv guillaumeKUnitigsAtLeast32bases_all.fasta.tmp guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta
rm -rf work1_mr
createSuperReadsForDirectory.perl -minreadsinsuperread 1 -l $KMER -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta -t $NUM_THREADS -mikedebug work1_mr superReadSequences.fasta.in 1> super1.err 2>&1
perl -ane 'push(@names,$F[0]);
END{
  open(FILE,"'work1_mr/superReadSequences.fasta.all'");
  while($line=<FILE>){
    if($line=~/^>/){
      chomp($line);
      print ">",$names[substr($line,1)],"\n";
    }else{
      print $line;
    }
  }
}' < work1_mr/superReadNames.txt > superReadSequences.named.fasta.tmp && mv superReadSequences.named.fasta.tmp superReadSequences.named.fasta
rm -f superReadSequences.fasta.in
touch .rerun
fi

SUPERREADS=superReadSequences.named.fasta
if [ ! -s superReadSequences.named.fasta ];then
error_exit "Error creating named super-reads file ";
fi

KUNITIGS=guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta
if [ ! -s $KUNITIGS ];then
error_exit "K-unitigs file $KUNITIGS not found!";
fi

REF_SPLIT="reference.split.fa"
if [ ! -s $REF_SPLIT ] || [ -e .rerun ];then
log "Preparing reference"
perl -ane '{
  if($F[0]=~/^>/){
    if(length($seq)>0){
      @f=split(/(N{1,})/,uc($seq)); 
      my $n=1;
      foreach $c(@f){
        if(not($c=~/^N/) && length($c)>0){
          $start=$n;
          $end=$n+length($c)-1;
          for(my $i=0;$i<length($c);$i+=10000000){
            if($i>0 && length($c)-$i<10000){
              print ">$rn:$start-$end:$i\n",substr($c,length($c)-10000,10000),"\n";
            }else{
              print ">$rn:$start-$end:$i\n",substr($c,$i,10000000),"\n"; 
            }
          }
        }
        $n+=length($c);
      }
    }
    $rn=substr($F[0],1);
    $seq="";
  }else{
    $seq.=$F[0];
  }
}END{
  @f=split(/(N{1,})/,uc($seq));
  my $n=1;
  foreach $c(@f){
    if(not($c=~/^N/) && length($c)>0){
      $start=$n;
      $end=$n+length($c)-1;
      for(my $i=0;$i<length($c);$i+=10000000){
        if($i>0 && length($c)-$i<10000){
          print ">$rn:$start-$end:$i\n",substr($c,length($c)-10000,10000),"\n";
        }else{
          print ">$rn:$start-$end:$i\n",substr($c,$i,10000000),"\n"; 
        }
      }
    }
    $n+=length($c);
  }
}'  $REF > $REF_SPLIT.tmp && mv $REF_SPLIT.tmp $REF_SPLIT
touch .rerun
fi

JF_SIZE=$(stat -c%s $KUNITIGS);

#the following two steps take the longest, so we run them in two parallel subshells

( if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
rm -f reconcile.success
log "Mega-reads pass 1"
if numactl --show 1> /dev/null 2>&1;then
numactl --interleave=all create_mega_reads -s $JF_SIZE -O 1.1 -e 5 -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 3000 -d $d  -r $SUPERREADS  -p $REF_SPLIT -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt || error_exit "Create_mega_reads failed"
else
create_mega_reads -s $JF_SIZE -O 1.1 -e 5 -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 3000 -d $d  -r $SUPERREADS  -p $REF_SPLIT -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt ||
 error_exit "Create_mega_reads failed"
fi
touch .rerun
fi

#this is different from joining the mega-reads, because we create scaffolds, not contigs; one scaffold per reference contig
if [ ! -s $COORDS.1.fa ] || [ -e .rerun ];then
log "Joining"
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt |join_mega_reads_trim.onepass.ref.pl 1>$COORDS.1.fa.tmp 2>/dev/null && mv $COORDS.1.fa.tmp $COORDS.1.fa
touch .rerun
fi ) &
PID1=$!

( if [ ! -s $COORDS.contigs.fa ] || [ -e .rerun ];then
rm -f reconcile.success
SR_FRG=$COORDS.sr.frg
log "Running preliminary assembly"
if [ ! -s $SR_FRG ];then
create_sr_frg.pl 65525 < $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta | fasta2frg.pl super > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG && rm -rf $CA
fi

rm -f $CA.log

echo "cnsConcurrency=$NUM_THREADS
cnsMinFrags=10000
unitigger=bogart
merylMemory=65536
ovlStoreMemory=65536
utgGraphErrorLimit=1000
utgMergeErrorLimit=1000
utgGraphErrorRate=0.015
utgMergeErrorRate=0.015
ovlCorrBatchSize=100000
ovlCorrConcurrency=8
frgCorrThreads=$NUM_THREADS
ovlThreads=2
ovlHashBlockLength=100000000
ovlRefBlockSize=1000000
ovlConcurrency=$NUM_THREADS
doFragmentCorrection=0
doOverlapBasedTrimming=0
doUnitigSplitting=0
doChimeraDetection=normal
merylThreads=$NUM_THREADS
doExtendClearRanges=0
cgwErrorRate=0.15
cgwMergeMissingThreshold=-1
cgwMergeFilterLevel=1
cgwDemoteRBP=0
cnsMaxCoverage=7
cnsReuseUnitigs=1" > runCA.spec

if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh $CA/5-consensus/consensus.sh
$CA_PATH/runCA -s runCA.spec -p genome -d $CA stopAfter=consensusAfterUnitigger $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
fi

if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
error_exit "Preliminary assembly failure, see $CA.log"
#else
#$CA_PATH/tigStore -g $CA/genome.gkpStore -t $CA/genome.tigStore 2 -U -d consensus -nreads 3 100000000 > $COORDS.unitigs.fa.tmp && mv $COORDS.unitigs.fa.tmp $COORDS.unitigs.fa
fi

if [ ! -e ${CA}/recompute_astat.success ];then
log "Recomputing A-stat"
recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
touch ${CA}/recompute_astat.success
fi

$CA_PATH/runCA -s runCA.spec -p genome -d $CA $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
cat $CA/9-terminator/genome.{ctg,deg}.fasta > $COORDS.contigs.fa.tmp && mv $COORDS.contigs.fa.tmp $COORDS.contigs.fa
touch .rerun
fi ) &
PID2=$!

wait $PID1 $PID2

#reconcile the reference "scaffolds" with the Illumina assembly
if [ ! -e reconcile.success ];then
rm -f final_assembly.success
log "Polishing reference contigs"
mkdir -p reconcile
(cd reconcile && polish_with_illumina_assembly.sh -r ../$COORDS.1.fa -q ../$COORDS.contigs.fa -t $NUM_THREADS -m 10000 1> /dev/null && \
splitScaffoldsAtNs.pl < $COORDS.1.fa.$COORDS.unitigs.fa.all.polished.deduplicated.fa > ../$COORDS.1.contigs.fa.tmp && mv ../$COORDS.1.contigs.fa.tmp ../$COORDS.1.contigs.fa ) && \
touch reconcile.success || error_exit "reconcile failed"
rm -f final_assembly.success
fi

#final assembly with Flye
if [ ! -e final_assembly.success ];then
log "Final assembly"
cat $COORDS.1.contigs.fa $COORDS.contigs.fa > $COORDS.subassemblies.fa && \
$MYPATH/flye -t $NUM_THREADS -i 0 --subassemblies $COORDS.subassemblies.fa  --kmer-size 25 -g $ESTIMATED_GENOME_SIZE -m 250 -o flye.$COORDS 1>flye.$COORDS.log 2>&1 && \
touch final_assembly.success || error_exit "Final assembly failure, see flye.$COORDS.log"
fi

if [ -e final_assembly.success ];then
log "Final output sequences are in flye.$COORDS/scaffolds.fasta"
ufasta n50 -a flye.$COORDS/scaffolds.fasta
fi
