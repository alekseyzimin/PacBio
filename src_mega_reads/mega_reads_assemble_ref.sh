#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
ESTIMATED_GENOME_SIZE=0
MER=17
B=25
d=0.05

#parsing arguments
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -m|--masurca_run_path)
    MASURCA_ASSEMBLY_WORK1_PATH="$2"
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
    -p|--pacbio)
    PACBIO="$2"
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
    echo "Usage: mega_reads_assemble.sh -m <path to MaSuRCA run work1 folder contents> -p <pacbio reads fasta> -a <path to the location of runCA in wgs-8.2 instalation>"
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

export PATH=$MYPATH:$CA_PATH:$PATH

if [ ! -e $PACBIO ];then
echo "PacBio reads file $PACBIO not found!";
exit 1 ;
fi

KUNITIGS=$MASURCA_ASSEMBLY_WORK1_PATH/../guillaumeKUnitigsAtLeast32bases_all.fasta
if [ ! -e $KUNITIGS ];then
echo "K-unitigs file $KUNITIGS not found!";
exit 1;
fi

KUNITIGLENGTHS=$MASURCA_ASSEMBLY_WORK1_PATH/kUnitigLengths.txt
if [ ! -e $KUNITIGLENGTHS ];then
echo "K-unitig lengths file $KUNITIGLENGTHS not found!";
exit 1;
fi

SUPERREADS=$MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta
if [ ! -e $SUPERREADS ];then
echo "Super reads file $SUPERREADS not found!";
exit 1;
else
if [ ! -e superReadSequences.named.fasta ];then
perl -ane 'push(@names,$F[0]);
	END{
	open(FILE,"'$SUPERREADS'");
	while($line=<FILE>){
		if($line=~/^>/){
			chomp($line);
			print ">",$names[substr($line,1)],"\n";
		}else{
			print $line;
	}
	}
}' < $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt > superReadSequences.named.fasta.tmp && mv superReadSequences.named.fasta.tmp superReadSequences.named.fasta
fi
if [ -s superReadSequences.named.fasta ];then
SUPERREADS=superReadSequences.named.fasta;
else
echo "Error creating named super-reads file from $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt and $SUPERREADS!";
rm superReadSequences.named.fasta
exit 1;
fi
fi

################setting parameters#########################
KMER=`perl -ane 'BEGIN{$min=10000}{if($F[1]<$min){$min=$F[1]}}END{print $min}' $KUNITIGLENGTHS`
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
JF_SIZE=`ls -l $SUPERREADS | perl -ane '{print $F[4]}'`
COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}

echo "Running mega-reads correction/assembly"
echo "Using mer size $MER for mapping, B=$B, d=$d"
echo "Using MaSuRCA files from $MASURCA_ASSEMBLY_WORK1_PATH, k-unitig mer $KMER"
echo "Using CA installation from $CA_PATH"
echo "Using $NUM_THREADS threads"
echo "Output prefix $COORDS"

rm -f .rerun

if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
echo "Mega-reads pass 1"
if numactl --show 1> /dev/null 2>&1;then
numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $PACBIO -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
else
create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $PACBIO -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
fi
touch .rerun
fi

if [ ! -s $COORDS.all.txt ] || [ -e .rerun ];then
echo "Refining alignments"
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt | split_matches_file.pl 1000 .matches && parallel "refine.sh $PACBIO $COORDS {1} $KMER" ::: .matches.* && cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && rm .matches.* && rm $COORDS.matches*.all.txt.tmp
#awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.mr.txt > $COORDS.all.txt.tmp && mv $COORDS.all.txt.tmp $COORDS.all.txt
touch .rerun
fi

if [ ! -s $COORDS.1.fa ] || [ -e .rerun ];then
echo "Joining"
join_mega_reads_trim.onepass.ref.pl < ${COORDS}.all.txt 1>$COORDS.1.fa.tmp 2>/dev/null && mv $COORDS.1.fa.tmp $COORDS.1.fa
touch .rerun
fi

if [ ! -s $COORDS.1.frg ] || [ -e .rerun ];then
echo "Generating assembly input files"
split_long_unitigs.pl refread < $COORDS.1.fa | make_mr_frg.pl mr 400 > $COORDS.1.frg.tmp && mv  $COORDS.1.frg.tmp  $COORDS.1.frg
rm -rf $CA
fi

TCOVERAGE=20
if [ $ESTIMATED_GENOME_SIZE -gt 1 ];then
MR_SIZE=$(stat -c%s "$COORDS.1.fa");
COVERAGE=$((MR_SIZE/ESTIMATED_GENOME_SIZE+4));
TCOVERAGE=`perl -e 'print int(int("'$COVERAGE'")/log(2)+1)'`
echo "Coverage threshold for splitting unitigs is $TCOVERAGE $BATOPTIONS"
SR_FRG=$COORDS.sr.frg
if [ ! -s $SR_FRG ];then
fasta2frg.pl sr 100 < $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG;
fi
fi

rm -f .rerun

echo "Running assembly"
runCA \
batOptions="-repeatdetect $TCOVERAGE $TCOVERAGE $TCOVERAGE" \
cnsConcurrency=$NUM_THREADS \
cnsMinFrags=1000 \
unitigger=bogart \
merylMemory=32768 \
ovlStoreMemory=32768 \
utgGraphErrorLimit=1000  \
utgMergeErrorLimit=1000 \
utgGraphErrorRate=0.025 \
utgMergeErrorRate=0.025 \
ovlCorrBatchSize=100000 \
ovlCorrConcurrency=4 \
frgCorrThreads=$NUM_THREADS \
mbtThreads=$NUM_THREADS \
ovlThreads=2 \
ovlMerThreshold=300 \
obtMerThreshold=400 \
ovlHashBlockLength=100000000 \
ovlRefBlockSize=1000000 \
ovlConcurrency=$NUM_THREADS \
doExtendClearRanges=0 \
doFragmentCorrection=0 \
doOverlapBasedTrimming=1 \
doUnitigSplitting=0 \
doChimeraDetection=normal \
-p genome -d $CA  \
merylThreads=$NUM_THREADS \
cnsReuseUnitigs=1 \
cgwMergeMissingThreshold=-1 \
cgwMergeFilterLevel=1 \
cgwDemoteRBP=0 \
cgwErrorRate=0.25 \
stopAfter=consensusAfterUnitigger \
$COORDS.1.frg $SR_FRG $OTHER_FRG 1> $CA.log 2>&1 && \
recompute_astat_superreads.sh genome $CA $PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  && \
runCA \
batOptions="-repeatdetect $TCOVERAGE $TCOVERAGE $TCOVERAGE" \
cnsConcurrency=$NUM_THREADS \
cnsMinFrags=1000 \
unitigger=bogart \
merylMemory=32768 \
ovlStoreMemory=32768 \
utgGraphErrorLimit=1000  \
utgMergeErrorLimit=1000 \
utgGraphErrorRate=0.025 \
utgMergeErrorRate=0.025 \
ovlCorrBatchSize=100000 \
ovlCorrConcurrency=4 \
frgCorrThreads=$NUM_THREADS \
mbtThreads=$NUM_THREADS \
ovlThreads=2 \
ovlMerThreshold=300 \
obtMerThreshold=400 \
ovlHashBlockLength=100000000 \
ovlRefBlockSize=1000000 \
ovlConcurrency=$NUM_THREADS \
doExtendClearRanges=0 \
doFragmentCorrection=0 \
doOverlapBasedTrimming=1 \
doUnitigSplitting=0 \
doChimeraDetection=normal \
-p genome -d $CA  \
merylThreads=$NUM_THREADS \
cnsReuseUnitigs=1 \
cgwMergeMissingThreshold=-1 \
cgwMergeFilterLevel=1 \
cgwDemoteRBP=0 \
cgwErrorRate=0.25 \
$COORDS.1.frg $SR_FRG $OTHER_FRG 1> $CA.log 2>&1 && \
echo "Assembly complete. Results are in $CA/9-terminator"

