#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -o pipefail
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export MYPATH
ESTIMATED_GENOME_SIZE=0
#minimum 50
MAX_GAP=1000
MER=15
B=17
d=0.029
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
PB_HC=30
KMER=0
#this is the batch size for grid execution
PBATCH_SIZE=2000000000
GRID_ENGINE="SGE"
QUEUE=""
USE_GRID=0
PACBIO=""
NANOPORE=""
NANOPORE_RNA=0
ONEPASS=0
OVLMIN_DEFAULT=250
POSTFIX=""
#MIN_PROPORTION="0.15"
#MIN_RADIUS="15"
MIN_PROPORTION="0.25"
MIN_RADIUS="400"
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
kill 0
exit 1
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
	-m|--masurca_run_path)
	    MASURCA_ASSEMBLY_WORK1_PATH="$2"
	    shift
	    ;;
	-t|--threads)
	    NUM_THREADS="$2"
	    shift
	    ;;
	-k|--kunitig_mer)
	    KMER="$2"
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
        -r|--rnaseq)
            NANOPORE="$2"
            NANOPORE_RNA=1
            shift
            ;;
        -n|--nanopore)
            NANOPORE="$2"
            shift
            ;;
        -C|--coverage)
            PB_HC="$2"
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
	-G|--use-grid)
	    USE_GRID="$2"
            shift
	    ;;
        -E|--grid-engine)
            GRID_ENGINE="$2"
            shift
            ;;
        -Pb|--pbatch-size)
            PBATCH_SIZE="$2"
            shift
            ;;
	-q|--queue)
	    QUEUE="$2"
	    shift
	    ;;
        -1|--onepass)
            ONEPASS=1;
            ;;
	-o|--other_frg)
	    OTHER_FRG="$2"
	    shift
	    ;;
        -F|--Flye)
            FLYE=1;
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
export PATH=$MYPATH:$CA_PATH:$PATH

if [ ! -s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta ];then
  error_exit "super reads file not found or size zero, you can try deleting $MASURCA_ASSEMBLY_WORK1_PATH folder and re-generating assemble.sh, also check if guillaumeKUnitigsAtLeast32bases_all.fasta is not empty";
fi

################setting parameters#########################
JF_SIZE=$(stat -L -c%s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta);
if [ $ESTIMATED_GENOME_SIZE -lt 1 ];then 
    error_exit "Estimated Genome Size is invalid or missing";
fi
COORDS=mr.$KMER.$MER.$B.$d

log "Running mega-reads correction/assembly"
log "Using mer size $MER for mapping, $KMER for k-unitigs, B=$B, d=$d"
log "Using $NUM_THREADS threads"
log "Output prefix $COORDS"

rm -f .rerun

if [ "$PACBIO" = "" ];then
    if [ "$NANOPORE" = "" ];then
	error_exit "must specify either PacBio or Nanopore input reads file with -p or -n"
    else
	LONGREADS=$NANOPORE
    fi
else
    LONGREADS=$PACBIO
fi

if [ ! -e $LONGREADS ];then
    error_exit "Long reads reads file $LONGREADS not found!";
fi

LONGREADS1="long_reads.fa";
if [ ! -s $LONGREADS1 ];then
  zcat -f $LONGREADS | $MYPATH/fastqToFasta.pl > $LONGREADS1.tmp && mv $LONGREADS1.tmp $LONGREADS1
fi

SUPERREADS=superReadSequences.named.fasta
KUNITIGS=$MASURCA_ASSEMBLY_WORK1_PATH/../guillaumeKUnitigsAtLeast32bases_all.fasta
if [ ! -s $SUPERREADS ];then
perl -ane 'push(@names,$F[0]);
END{
  open(FILE,"'$MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta.all'");
  while($line=<FILE>){
    if($line=~/^>/){
      chomp($line);
      print ">",$names[substr($line,1)],"\n";
    }else{
      print $line;
    }
  }
}' < $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt > $SUPERREADS.tmp && mv $SUPERREADS.tmp $SUPERREADS || error_exit "failed to create named super-reads file";
fi

if [ ! -s $SUPERREADS ];then
    error_exit "Error creating named super-reads file ";
fi

if [ ! -s $KUNITIGS ];then
    error_exit "K-unitigs file $KUNITIGS not found; failed to create named super-reads file";
fi

KUNITIGLENGTHS=$MASURCA_ASSEMBLY_WORK1_PATH/kUnitigLengths.txt
if [ ! -s $KUNITIGLENGTHS ];then
    error_exit "K-unitig lengths file $KUNITIGLENGTHS not found!";
fi


if numactl --show 1> /dev/null 2>&1;then #if numactl then interleave memory
  numactl --interleave=all create_mega_reads -L 50 -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 10000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed"
else
  create_mega_reads -L 50 -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 10000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed"
fi

if [ ! -s ${COORDS}.all.txt ] || [ -e .rerun ];then
    log "Refining alignments"
    NUM_LONGREADS_READS_PER_BATCH=`grep --text '^>'  $LONGREADS1 | wc -l | awk '{bs=int($1/1024);if(bs<1000){bs=1000};if(bs>100000){bs=100000};}END{print bs}'` 
    awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt | add_pb_seq.pl $LONGREADS1 | split_matches_file.pl $NUM_LONGREADS_READS_PER_BATCH .matches && ls .matches.* | xargs -P $NUM_THREADS -I % refine.sh $COORDS % $KMER && cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && rm .matches.* && rm $COORDS.matches*.all.txt.tmp 
    touch .rerun
fi

#output -- just the corrected chunks
grep -v '^>'  $COORDS.all.txt |awk '{print ">"$6":"$1"-"$2"\n"$7}' > $COORDS.transcripts.fa
log "OUtput corrected transcripts are in $COORDS.transcripts.fa"

