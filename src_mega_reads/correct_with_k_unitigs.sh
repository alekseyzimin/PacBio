#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
export NUM_THREADS=4
export PKMER=19
export LONGREADS=""
export LONGREADS1=""
export PKUNITIGS=k_unitigs.l20.fa
export COVERAGE=1000
export COORDS="pCorrect"
export ESTIMATED_GENOME_SIZE=0
export LONGREADS="in.fa"
export LONGREADS1="out.fa"
export ILLUMINA="illumina.fa"
USAGE="Usage:  correct_with_k_unitigs.sh -k <k-unitig k-mer, default:19> -e <estimated genome size, MANDATORY> -t <number of threads, default:4> -l <input long reads file name, default:in.fa> -i <illumina reads file name, default illumina.fa> -o <output long reads file name, default:out.fa>"
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
kill -9 0
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
if [[ $# -eq 0 ]];then
echo $USAGE
exit 0
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -k|--kmer)
            export PKMER="$2"
            shift
            ;;
        -t|--threads)
            export NUM_THREADS="$2"
            shift
            ;;
        -e)
            export ESTIMATED_GENOME_SIZE="$2"
            shift
            ;;
        -l|--long_reads)
            export LONGREADS="$2"
            shift
            ;;
        -o|--output)
            export LONGREADS1="$2";
            shift
            ;;
        -i|--illumina)
            export ILLUMINA="$2";
            shift
            ;;
        -c|--coverage)
            export COVERAGE="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo $USAGE
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

#check arguments
if [ $ESTIMATED_GENOME_SIZE -lt 1 ];then
  echo $USAGE
  error_exit "Must specify estimated genome size that is greater than zero"
fi

if [ ! -s $LONGREADS ];then
  echo $USAGE
  error_exit "Input file of long reads $LONGREADS does not exist or size zero"
fi

if [ ! -s $ILLUMINA ];then
  echo $USAGE
  error_exit "Input file of illumina reads $ILLUMINA does not exist or size zero"
fi

export PKUNITIGS=k_unitigs.l$(($PKMER+1)).fa

if [ ! -s $PKUNITIGS ];then
  log "Creating k-unitigs for k=$PKMER"
  $MYPATH/create_k_unitigs_large_k2 -c $(($PKMER-1)) -t $NUM_THREADS -m $PKMER -n $(($ESTIMATED_GENOME_SIZE*2)) -l $(($PKMER+1)) $ILLUMINA  | \
  grep --text -v '^>' | \
  perl -ane '{$seq=$F[0]; $F[0]=~tr/ACTGactg/TGACtgacc/;$revseq=reverse($F[0]); $h{($seq ge $revseq)?$seq:$revseq}=1;}END{$n=0;foreach $k(keys %h){print ">",$n++," length:",length($k),"\n$k\n"}}' > $PKUNITIGS.tmp && \
  mv $PKUNITIGS.tmp $PKUNITIGS
fi

if [ ! -s $COORDS.pcorrected.fa ];then
  log "Pre-correcting long reads"
  rm -f $COORDS.pcorrected.*.fa
  echo "#!/bin/bash" > correct_with_k_unitigs.sh
  echo "$MYPATH/create_mega_reads -s \$(($ESTIMATED_GENOME_SIZE*2)) -m $PKMER --psa-min 12 --stretch-cap 10000 -k $PKMER -u $PKUNITIGS -t $NUM_THREADS -B 1 --max-count 5000 -d 0.01 -r <(awk '{if(\$1 ~ /^>/) print \$1\"F\"; else print \$1}'  $PKUNITIGS)  -p <(zcat -f $LONGREADS | $MYPATH/fastqToFasta.pl | awk 'BEGIN{rn=0;}{if(\$1 ~ /^>/){print \">\"rn;rn++;}else{print \$1}}') -L $PKMER -o /dev/stdout 2>create_mega-reads.err |\\" >> correct_with_k_unitigs.sh
  echo "$MYPATH/add_pb_seq.pl <(zcat -f $LONGREADS | $MYPATH/fastqToFasta.pl | awk 'BEGIN{rn=0;}{if(\$1 ~ /^>/){print \">\"rn;rn++;}else{print \$1}}') |\\" >> correct_with_k_unitigs.sh
  echo "$MYPATH/ufasta split -i /dev/stdin \\" >> correct_with_k_unitigs.sh
  for i in $(seq 1 $(($NUM_THREADS/16+2)));do
    echo ">($MYPATH/correct_with_k_unitigs_fast.pl $PKMER 0.0 1>$COORDS.pcorrected.$i.fa 2>/dev/null) \\" >> correct_with_k_unitigs.sh
  done
  echo "&& cat $COORDS.pcorrected.*.fa | awk '{if(\$1!~/^>/ && \$1~/>/){split(\$1,a,\">\");print a[1];if(a[2]!=\"\") print \">\"a[2];}else{print \$0}}'  > $COORDS.pcorrected.fa.tmp && mv $COORDS.pcorrected.fa.tmp $COORDS.pcorrected.fa && rm -f $COORDS.pcorrected.*.fa" >> correct_with_k_unitigs.sh
  chmod 0755 correct_with_k_unitigs.sh && ./correct_with_k_unitigs.sh
fi

if [ -s $COORDS.pcorrected.fa ];then
  $MYPATH/ufasta extract -f <(paste <($MYPATH/ufasta sizes -H  $COORDS.pcorrected.fa) <(tr -d acgt < $COORDS.pcorrected.fa |$MYPATH/ufasta sizes -H ) | \
  sort -nrk4 -S 10% | awk '{n+=$2;if(n<int('$ESTIMATED_GENOME_SIZE')*int('$COVERAGE')) print $1}' ) $COORDS.pcorrected.fa > $LONGREADS1.tmp && \
  mv $LONGREADS1.tmp $LONGREADS1 && \
  rm -f $COORDS.pcorrected.fa correct_with_k_unitigs.sh $PKUNITIGS
else
  error_exit "Failed to pre-correct long reads, please check your input files $LONGREADS and $ILLUMINA and the other parameters"
fi
log "Pre-corrected reads are in $LONGREADS1"
