#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
FLYE_PATH="$MYPATH/../Flye/bin";
CA_PATH="$MYPATH/../CA8/Linux-amd64/bin";
ESTIMATED_GENOME_SIZE=0
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
MER=17
B=13
d=0.03
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

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0 
exit 1
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
error_exit "runCA not found at $CA_PATH, please check your MaSuRCA install";
fi

if [ ! -e $FLYE_PATH/flye ];then
error_exit "flye not found at $FLYE_PATH, please check your MaSuRCA install";
fi

export PATH=$CA_PATH:$MYPATH:$PATH:$FLYE_PATH

if [ ! -e $REF ];then
error_exit "Reference sequence file $REF not found!";
fi

COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}

log "Running synteny-guided assembly"
log "Using mer size $MER for mapping, B=$B, d=$d"
log "Using CA installation from $CA_PATH"
log "Using $NUM_THREADS threads"
log "Output prefix $COORDS"

SUPERREADS=superReadSequences.named.fasta
KUNITIGS=guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta

#first we re-create k-unitigs and super reads with smaller K
if [ ! -e reduce-super-reads.success ];then
rm -f mega-reads.success
log "Reducing super-read k-mer size"
awk 'BEGIN{n=0}{if($1~/^>/){}else{print ">sr"n"\n"$0;n+=2;}}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta > superReadSequences.fasta.in && \
create_k_unitigs_large_k -q 1 -c $(($KMER-1)) -t $NUM_THREADS -m $KMER -n $(($ESTIMATED_GENOME_SIZE*2)) -l $KMER -f `perl -e 'print 1/'$KMER'/1e5'` $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta.all  | grep --text -v '^>' | perl -ane '{$seq=$F[0]; $F[0]=~tr/ACTGactg/TGACtgac/;$revseq=reverse($F[0]); $h{($seq ge $revseq)?$seq:$revseq}=1;}END{$n=0;foreach $k(keys %h){print ">",$n++," length:",length($k),"\n$k\n"}}' > $KUNITIGS.tmp && mv $KUNITIGS.tmp $KUNITIGS && \
rm -rf work1_mr && \
createSuperReadsForDirectory.perl -minreadsinsuperread 1 -l $KMER -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile $KUNITIGS -t $NUM_THREADS -mikedebug work1_mr superReadSequences.fasta.in 1> super1.err 2>&1 && \
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
}' < work1_mr/superReadNames.txt > $SUPERREADS.tmp && mv $SUPERREADS.tmp $SUPERREADS && touch reduce-super-reads.success && rm -f superReadSequences.fasta.in
if [ ! -e reduce-super-reads.success ];then
  error_exit "Failed to reduce super-reads k-mer size, check super1.err"
fi
fi

if [ ! -s $SUPERREADS ];then
error_exit "Error creating named super-reads file";
fi

if [ ! -s $KUNITIGS ];then
error_exit "K-unitigs file $KUNITIGS not found, failed to reduce super-reads k-mer size";
fi

REF_SPLIT="reference.split.fa"
if [ ! -e prepare.success ];then
rm -f mega-reads.success
log "Preparing reference"
perl -ane 'BEGIN{$chunk_size=25000000;}{
  if($F[0]=~/^>/){
    if(length($seq)>0){
      @f=split(/(N{1,})/,uc($seq)); 
      my $n=1;
      foreach $c(@f){
        if(not($c=~/^N/) && length($c)>0){
          $start=$n;
          $end=$n+length($c)-1;
          for(my $i=0;$i<length($c);$i+=$chunk_size){
            if($i>0 && length($c)-$i<10000){
              print ">$rn:$start-$end:$i\n",substr($c,length($c)-10000,10000),"\n";
            }else{
              print ">$rn:$start-$end:$i\n",substr($c,$i,$chunk_size),"\n";
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
          print ">$rn:$start-$end:$i\n",substr($c,$i,$chunk_size),"\n"; 
        }
      }
    }
    $n+=length($c);
  }
}'  $REF > $REF_SPLIT.tmp && mv $REF_SPLIT.tmp $REF_SPLIT &&  touch prepare.success
fi

JF_SIZE=$(stat -c%s $KUNITIGS);

#the following two steps take the longest, so we run them in two parallel subshells

( if [ ! -e mega-reads.success ];then
rm -f reconcile.success
rm -f join.success
log "Building synteny-reads"
if numactl --show 1> /dev/null 2>&1;then
numactl --interleave=all create_mega_reads -s $JF_SIZE -O 1.1 -e 5 -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 3000 -d $d  -r $SUPERREADS  -p $REF_SPLIT -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt && touch mega-reads.success || error_exit "Create_mega_reads failed"
else
create_mega_reads -s $JF_SIZE -O 1.1 -e 5 -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 3000 -d $d  -r $SUPERREADS  -p $REF_SPLIT -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt && touch mega-reads.success || error_exit "Create_mega_reads failed"
fi
fi

#this is different from joining the mega-reads, because we create scaffolds, not contigs; one scaffold per reference contig
if [ ! -e join.success ];then
log "Joining"
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt |join_mega_reads_trim.onepass.ref.pl 1>$COORDS.1.fa.tmp 2>/dev/null && mv $COORDS.1.fa.tmp $COORDS.1.fa && touch join.success
fi ) &
PID1=$!

( if [ ! -e CA.success ];then
rm -f reconcile.success
SR_FRG=superReads.frg
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
fi

if [ ! -e ${CA}/recompute_astat.success ];then
log "Recomputing A-stat"
recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
touch ${CA}/recompute_astat.success
fi

$CA_PATH/runCA -s runCA.spec -p genome -d $CA $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
ln -s $CA/9-terminator/genome.ctg.fasta CA.contigs.fa && touch CA.success
fi ) &
PID2=$!

wait $PID1 $PID2

#final assembly with Flye
if [ ! -e final_assembly.success ];then
log "Final assembly"
rm -f polish.success
cat <(splitScaffoldsAtNs.pl < $COORDS.1.fa) CA.contigs.fa > $COORDS.subassemblies.fa.tmp && mv $COORDS.subassemblies.fa.tmp $COORDS.subassemblies.fa && \
$FLYE_PATH/flye -t $NUM_THREADS -i 0 --subassemblies $COORDS.subassemblies.fa  --kmer-size 25 -g $ESTIMATED_GENOME_SIZE -m 250 -o flye.$COORDS 1>flye.$COORDS.log 2>&1 && \
touch final_assembly.success || error_exit "Final assembly failure, see flye.$COORDS.log"
fi

if [ -e final_assembly.success ];then
ufasta extract -f <(ufasta sizes -H flye.$COORDS/assembly.fasta | awk '{if($2>=1000) print $1}') flye.$COORDS/assembly.fasta > final_assembly.fasta.tmp && \
mv final_assembly.fasta.tmp final_assembly.fasta && \
ufasta n50 -a final_assembly.fasta
fi

#if [ ! -e polish.success ];then
#log "Polishing assembly"
#(mkdir -p polish.$COORDS && \
#cd polish.$COORDS && \
#ufasta extract -f <(ufasta sizes -H ../flye.$COORDS/assembly.fasta | awk '{if($2>=1000) print $1}') ../flye.$COORDS/assembly.fasta > filtered_assembly.fasta.tmp && mv filtered_assembly.fasta.tmp filtered_assembly.fasta && \
#polca.sh -a filtered_assembly.fasta -r ../pe.cor.fa -t $NUM_THREADS -m 1G && \
#cd ..) && touch polish.success || error_exit "Polishing assembly failure, see files in polish.$COORDS/"
#fi

#if [ -e polish.success ];then
#ln -sf polish.$COORDS/filtered_assembly.fasta.PolcaCorrected.fa final_assembly.fasta && \
#log "Success! Final output sequences are in final_assembly.fasta" && \
#ufasta n50 -a final_assembly.fasta
#fi



