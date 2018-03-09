#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
ESTIMATED_GENOME_SIZE=0
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`

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

export PATH=$CA_PATH:$MYPATH:$PATH

if [ ! -e $REF ];then
echo "Reference reads file $REF not found!";
exit 1 ;
fi
################setting parameters#########################
MER=17
B=20
d=0.05
KMER=41
COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}

echo "Running mega-reads correction/assembly"
echo "Using mer size $MER for mapping, B=$B, d=$d"
echo "Using MaSuRCA files from $MASURCA_ASSEMBLY_WORK1_PATH, k-unitig mer $KMER"
echo "Using CA installation from $CA_PATH"
echo "Using $NUM_THREADS threads"
echo "Output prefix $COORDS"

rm -f .rerun

#first we re-create k-unitigs and super reads with smaller K
if [ ! -e superReadSequences.named.fasta ];then
echo "Reducing super-read k-mer size"
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
fi

SUPERREADS=superReadSequences.named.fasta
if [ ! -s superReadSequences.named.fasta ];then
echo "Error creating named super-reads file ";
exit 1;
fi

KUNITIGS=guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta
if [ ! -s $KUNITIGS ];then
echo "K-unitigs file $KUNITIGS not found!";
exit 1;
fi

REF_SPLIT="reference.split.fa"
if [ ! -s $REF_SPLIT ] || [ -e .rerun ];then
echo "Preparing reference"
perl -ane '{
  if($F[0]=~/^>/){
    if(length($seq)>0){
      @f=split(/(N{1,})/,uc($seq)); 
      my $n=1;
      foreach $c(@f){
        if(not($c=~/^N/) && length($c)>0){
          $start=$n;
          $end=$n+length($c)-1;
          print ">$rn:$start-$end\n$c\n";
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
      print ">$rn:$start-$end\n$c\n";
    }
    $n+=length($c);
  }
}'  $REF | create_sr_frg.pl 10000000 ref  > $REF_SPLIT.tmp && mv $REF_SPLIT.tmp $REF_SPLIT
fi

JF_SIZE=$(stat -c%s $KUNITIGS);

if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
echo "Mega-reads pass 1"
if numactl --show 1> /dev/null 2>&1;then
numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 3000 -d $d  -r $SUPERREADS  -p $REF_SPLIT -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
else
create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 3000 -d $d  -r $SUPERREADS  -p $REF_SPLIT -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
fi
touch .rerun
fi

if [ ! -s $COORDS.all.txt ] || [ -e .rerun ];then
echo "Refining alignments"
NUM_LONGREADS_READS_PER_BATCH=`grep --text '^>'  $REF_SPLIT | wc -l | awk '{bs=int($1/100);if(bs<100){bs=100};if(bs>100000){bs=100000};}END{print bs}'`
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt | add_pb_seq.pl $REF_SPLIT | split_matches_file.pl $NUM_LONGREADS_READS_PER_BATCH .matches && ls .matches.* | xargs -P $NUM_THREADS -I % refine.sh $COORDS % $KMER && cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && rm .matches.* && rm $COORDS.matches*.all.txt.tmp
#awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt | add_pb_seq.pl $REF_SPLIT > .matches.0 && refine.sh $COORDS .matches.0 $KMER && mv $COORDS.matches.0.all.txt.tmp $COORDS.all.txt && rm .matches.0
touch .rerun
fi

#this is different from joining the mega-reads, because we create scaffolds, not contigs; one scaffold per reference contig
if [ ! -s $COORDS.1.fa ] || [ -e .rerun ];then
echo "Joining"
join_mega_reads_trim.onepass.ref.pl < ${COORDS}.all.txt 1>$COORDS.1.fa.tmp 2>/dev/null && mv $COORDS.1.fa.tmp $COORDS.1.fa
touch .rerun
fi

#now we attempt to close gaps in reference assisted scaffolds
if [ ! -s $COORDS.1.gapclose.fa ] || [ -e .rerun ];then
echo "Closing gaps in reference assisted scaffolds"
(mkdir -p ref_gapclose && \
cd ref_gapclose && \
closeGapsInScaffFastaFile.perl --split 1 --max-reads-in-memory 1000000000 -s $(($ESTIMATED_GENOME_SIZE*5)) --scaffold-fasta-file ../$COORDS.1.fa  --reads-file ../pe.cor.fa --output-directory gapclose.tmp --min-kmer-len 19 --max-kmer-len 127 --num-threads $NUM_THREADS --contig-length-for-joining 300 --contig-length-for-fishing 450 --reduce-read-set-kmer-size 25 1>gapClose.err 2>&1 && \
mv gapclose.tmp/genome.scf.fasta ../$COORDS.1.gapclose.fa ) || error_exit "reference gapclose failed"
fi

SR_FRG=$COORDS.sr.frg
TCOVERAGE=20
if [ $ESTIMATED_GENOME_SIZE -gt 1 ];then
echo "Generating assembly input files"
if [ ! -s $SR_FRG ];then
create_sr_frg.pl 65525 < $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta | fasta2frg.pl super > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG && rm -rf $CA
fi
TCOVERAGE=`ls $SR_FRG $OTHER_FRG 2>/dev/null | xargs stat -c%s | awk '{n+=$1}END{print int(n/int('$ESTIMATED_GENOME_SIZE')/1.7+1)}'`;
fi

rm -f .rerun
rm -f $CA.log

OVLMIN=`head -n 100000 $SR_FRG $OTHER_FRG 2>/dev/null | grep -A 1 '^seq:' |grep -v '^seq:' | grep -v '\-\-' | awk 'BEGIN{minlen=100000}{if(length($1)<minlen && length($1)>=64) minlen=length($1);}END{if(minlen>=250) print "250"; else print minlen-5;}'`

batOptions="-repeatdetect $TCOVERAGE $TCOVERAGE $TCOVERAGE -el $OVLMIN "

echo "Coverage threshold for splitting unitigs is $TCOVERAGE minimum ovl $OVLMIN"

echo "batOptions=$batOptions
cnsConcurrency=$NUM_THREADS
cnsMinFrags=10000
unitigger=bogart
merylMemory=65536
ovlStoreMemory=65536
utgGraphErrorLimit=1000
utgMergeErrorLimit=1000
utgGraphErrorRate=0.02
utgMergeErrorRate=0.02
ovlCorrBatchSize=100000
ovlCorrConcurrency=4
frgCorrThreads=$NUM_THREADS
mbtThreads=$NUM_THREADS
ovlThreads=2
ovlHashBlockLength=100000000
ovlRefBlockSize=1000000
ovlConcurrency=$NUM_THREADS
doFragmentCorrection=1
doOverlapBasedTrimming=0
doUnitigSplitting=0
doChimeraDetection=normal
merylThreads=$NUM_THREADS
doExtendClearRanges=0
cgwErrorRate=0.15
cgwMergeMissingThreshold=-1
cgwMergeFilterLevel=1
cgwDemoteRBP=0
cnsReuseUnitigs=1" > runCA.spec

echo "Running assembly"
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
#need to start from the beginning
$CA_PATH/runCA -s runCA.spec -p genome -d $CA stopAfter=consensusAfterUnitigger $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
fi

#at athis point we assume that the unitig consensus is done
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
echo "CA failure, see $CA.log"
exit;
fi

if [ ! -e ${CA}/recompute_astat.success ];then
echo "Recomputing A-stat"
recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
touch ${CA}/recompute_astat.success
fi

#we start from here if the scaffolder has been run or continue here  
$CA_PATH/runCA -s runCA.spec -p genome -d $CA $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1

#reconcile the reference "scaffolds" with the Illumina assembly
if [ ! -e reconcile.success ];then
echo "Polishing reference contigs"
mkdir -p reconcile
(cd reconcile && polish_with_illumina_assembly.sh -r ../$COORDS.1.gapclose.fa -q ../$CA/9-terminator/genome.ctg.fasta -t $NUM_THREADS 1> /dev/null && perl -ane '{
  if($F[0]=~/^>/){
    if(length($seq)>0){
      @f=split(/(N{1,})/,uc($seq)); 
      my $n=1;
      foreach $c(@f){
        if(not($c=~/^N/) && length($c)>0){
          $start=$n;
          $end=$n+length($c)-1;
          print ">$rn:$start-$end\n$c\n";
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
      print ">$rn:$start-$end\n$c\n";
    }
    $n+=length($c);
  }
}'  $COORDS.1.gapclose.fa.genome.ctg.fasta.all.polished.deduplicated.fa > ../$COORDS.1.contigs.fa) && touch reconcile.success || error_exit "reconcile failed"
rm -f merge.success
fi

#now we merge
if [ ! -e merge.success ];then
echo "Merging reference contigs"
mkdir -p final_merge && \
(cd final_merge && merge_contigs.sh -r ../$COORDS.1.contigs.fa -q ../$CA/9-terminator/genome.ctg.fasta -t $NUM_THREADS 1>/dev/null && mv  $COORDS.1.contigs.fa.genome.ctg.fasta.merged.fa ../contigs.final.fasta) && touch merge.success || error_exit "final merge failed"
fi

echo "Final output contigs are in contigs.final.fasta"
ufasta n50 -A -S -N50 -C contigs.final.fasta
