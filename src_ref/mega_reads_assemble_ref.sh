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
MER=17
B=25
d=0.05
KMER=`awk 'BEGIN{min=10000}{if($2<min) min=$2}END{print min}' $KUNITIGLENGTHS`
JF_SIZE=$(stat -c%s $SUPERREADS);
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
COVERAGE=$((MR_SIZE/ESTIMATED_GENOME_SIZE+1));
if [ $COVERAGE -le 5 ];then
echo "Coverage of the mega-reads less than 5 -- using the super reads as well";
SR_FRG=$COORDS.sr.frg
if [ ! -s $SR_FRG ];then
awk '{if($0 ~ /^>/) print $0":super-read"; else print $0}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta | fasta2frg.pl sr 200 > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG;
fi
fi
COVERAGE=`ls $SR_FRG $COORDS.1.frg $COORDS.1.mates.frg $OTHER_FRG 2>/dev/null | xargs stat -c%s | awk '{n+=$1}END{print int(n/int('$ESTIMATED_GENOME_SIZE')/1.7+1)}'`;
TCOVERAGE=$COVERAGE;
fi

rm -f .rerun
rm -f $CA.log

OVLMIN=`head -n 100000 $SR_FRG $COORDS.1.frg $COORDS.1.mates.frg $OTHER_FRG 2>/dev/null | grep -A 1 '^seq:' |grep -v '^seq:' | grep -v '\-\-' | awk 'BEGIN{minlen=100000}{if(length($1)<minlen && length($1)>=64) minlen=length($1);}END{if(minlen>=250) print "250"; else print minlen-5;}'`

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
doOverlapBasedTrimming=1
doUnitigSplitting=0
doChimeraDetection=normal
merylThreads=$NUM_THREADS
doExtendClearRanges=1
cgwErrorRate=0.15
cgwMergeMissingThreshold=-1
cgwMergeFilterLevel=1
cgwDemoteRBP=0
cnsReuseUnitigs=1" > runCA.spec

echo "Running assembly"
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
#need to start from the beginning
runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
runCA -s runCA.spec -p genome -d $CA  stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
fi

#at athis point we assume that the unitig consensus is done
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
echo "Unitig consensus failure"
exit;
fi

#recompute astat if low pacbio coverage
if [ $MCOVERAGE -le 5 ]; then
if [ ! -e ${CA}/recompute_astat.success ];then
recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
fi
fi

#we start from here if the scaffolder has been run or continue here  
runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopAfter=consensusAfterScaffolder $COORDS.1.frg $SR_FRG $COORDS.1.mates.frg $OTHER_FRG 1>> $CA.log 2>&1
rm -rf $CA/8-consensus/*.success $CA/8-consensus/consensus.sh
runCA -s runCA.spec -p genome -d $CA  stopAfter=consensusAfterScaffolder $COORDS.1.frg $SR_FRG $COORDS.1.mates.frg $OTHER_FRG 1>> $CA.log 2>&1 && echo "Assembly complete. Results are in $CA/9-terminator"

