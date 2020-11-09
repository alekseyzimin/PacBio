#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -o pipefail
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export MYPATH
MER=15
B=15
d=0.03
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
ASSEMBLY=""
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
	-M|--alignment_mer)
	    MER="$2"
	    shift
	    ;;
        -l|--l_mer)
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
	-A|--Assembly)
	    ASSEMBLY="$2"
	    shift
	    ;;
	-v|--verbose)
	    set -x
	    ;;    
	-h|--help|-u|--usage)
	    echo "Usage: mega_reads_assemble_polish.sh -m <path to MaSuRCA run work1 folder contents> -A <fasta assembly to be polshed> -t <number of threads>"
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
export PATH=$MYPATH:$PATH

if [ ! -s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta ];then
  error_exit "super reads file not found or size zero, you can try deleting $MASURCA_ASSEMBLY_WORK1_PATH folder and re-generating assemble.sh, also check if guillaumeKUnitigsAtLeast32bases_all.fasta is not empty";
fi

################setting parameters#########################
JF_SIZE=$(stat -L -c%s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta);
COORDS=mr.$KMER.$MER.$B.$d
log "Running mega-reads polishing on $ASSEMBLY"
log "Using mer size $MER for mapping, B=$B, d=$d"
log "Estimated Ploidy $PLOIDY"
log "Using $NUM_THREADS threads"
log "Output prefix $COORDS"
log "Using kmer size $KMER for k-unitigs"
rm -f .rerun
###############removing redundant subreads or reducing the coverage by picking the longest reads##############################
if [ ! -s $ASSEMBLY ];then
        error_exit "$ASSEMBLY file is missing or size zero!"
fi

SUPERREADS=superReadSequences.named.fasta
KUNITIGS=$MASURCA_ASSEMBLY_WORK1_PATH/../guillaumeKUnitigsAtLeast32bases_all.fasta

if [ ! -e $COORDS.namedSuperReads.polish.success ];then
perl -ane 'push(@names,$F[0]);
END{
  open(FILE,"'$MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta'");
  while($line=<FILE>){
    if($line=~/^>/){
      chomp($line);
      print ">",$names[substr($line,1)],"\n";
    }else{
      print $line;
    }
  }
}' < $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt > superReadSequences.named.fasta.tmp && mv superReadSequences.named.fasta.tmp superReadSequences.named.fasta && touch $COORDS.namedSuperReads.polish.success  || error_exit "failed to create named super-reads file";
rm -f $COORDS.mega-reads.polish.success
fi

if [ ! -s $KUNITIGS ];then
    error_exit "K-unitigs file $KUNITIGS empty or not found!";
fi

if [ ! -e $COORDS.mega-reads.polish.success ];then
log "Creating mega-reads"
$MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $ASSEMBLY -o $COORDS.txt.tmp 2>create_mega_reads.err && mv $COORDS.txt.tmp $COORDS.txt && touch $COORDS.mega-reads.polish.success  || error_exit "mega-reads failed, see create_mega_reads.err"
rm -f $COORDS.polish.polish.success
fi

if [ ! -e $COORDS.polish.polish.success ];then
log "Realignment and polishing"
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt | add_pb_seq.pl $ASSEMBLY | refine_alignments.pl $COORDS > $COORDS.delta.tmp && mv $COORDS.delta.tmp $COORDS.delta && \
show-coords -lcHr  $COORDS.delta| reconcile_consensus.pl $ASSEMBLY t.$COORDS.maximal_mr.fa > $COORDS.polished.fa.tmp && mv $COORDS.polished.fa.tmp $COORDS.polished.fa && touch $COORDS.polish.polish.success || error_exit "polishing failed, please check input files and available disk space"
fi

if [ -e $COORDS.polish.polish.success ];then
log "Success! Polished assembly is in $COORDS.polished.fa"
ufasta n50 -a $COORDS.polished.fa
fi

