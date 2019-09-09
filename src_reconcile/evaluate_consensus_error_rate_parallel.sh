#!/bin/bash

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
set -e
NUM_THREADS=4
MEM=16000000000
FIX=0

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            export NUM_THREADS="$2"
            shift
            ;;
        -f|--fix)
            FIX=1
            ;;
        -a|--assembly)
            export ASM="$2"
            shift
            ;;
        -r|--reads)
            READS="$2";
            shift
            ;;
        -m|--memory)
            export MEM="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage:  evaluate_consensus_error_rate.sh -a <assembly contigs or scaffolds> -r <Illumina reads fastq> -t <number of threads [-f] optional:fix errors that are found>"
            echo "Must have bwa, samtools and freebayes available on the PATH"
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [ ! -e $ASM ] || [ ! -e $READS ];then
echo "assembly file $ASM or Illumina reads file $READS not found!"
exit 1
fi

export BASM=`basename $ASM`
export BWA=`which bwa`
export FREEBAYES=`which freebayes`
export SAMTOOLS=`which samtools`
rm -f bwa.err samtools.err

if [ ! -e $BASM.index.success ];then 
echo "Creating BWA index for $ASM"
rm -f $BASM.map.success
$BWA index $ASM -p $BASM.bwa 2>>bwa.err && touch $BASM.index.success
fi

if [ ! -e $BASM.map.success ];then
echo "Aligning reads to $ASM"
rm -f $BASM.sort.success
$BWA mem -SP -t $NUM_THREADS $BASM.bwa $READS 2>>bwa.err |samtools view -bhS /dev/stdin 1>$BASM.unSorted.bam 2>>samtools.err && touch $BASM.map.success
fi

if [ ! -e $BASM.sort.success ];then
echo "Sorting and indexing alignment file"
rm -f $BASM.vc.success
$SAMTOOLS sort -m $MEM  $BASM.unSorted.bam $BASM.alignSorted 2>>samtools.err && \
$SAMTOOLS index $BASM.alignSorted.bam 2>>samtools.err && touch  $BASM.sort.success
fi

#here we are doing variant calling in parallel, per input contig/scaffold
if [ ! -e $BASM.vc.success ];then
rm -f  $BASM.fix.success
echo "Calling variants"
grep '^>' $ASM |awk '{print substr($1,2)}' > $BASM.names
mkdir -p $BASM.work
#subshell
(
  cd $BASM.work
  CONTIGS=`wc -l ../$BASM.names | awk '{print $1}'`;
  BATCH_SIZE=$(($CONTIGS / $NUM_THREADS+1));
  if [ $BATCH_SIZE -gt 1000 ];then
    BATCH_SIZE=1000
  fi
  BATCH=1;
  echo "Processing $BATCH_SIZE scaffold(s) per batch"
  LIST="";
  INDEX=1;
  for f in $(cat ../$BASM.names);do
    LIST="$LIST $f"
    let INDEX=$INDEX+1;
      if [ $INDEX -gt $BATCH_SIZE ];then
        if [ ! -e $BATCH.extract.success ];then
          $SAMTOOLS view -h ../$BASM.alignSorted.bam $LIST 2>>samtools.err |$SAMTOOLS view -S -b /dev/stdin 2>>samtools.err 1> $BATCH.alignSorted.bam && touch $BATCH.extract.success
        fi
          LIST=""
          INDEX=1
          let BATCH=$BATCH+1
      fi
  done
  if [ $INDEX -gt 1 ];then
    if [ ! -e $BATCH.extract.success ];then
      $SAMTOOLS view -h ../$BASM.alignSorted.bam $LIST 2>>samtools.err |$SAMTOOLS view -S -b /dev/stdin 2>>samtools.err 1> $BATCH.alignSorted.bam && touch $BATCH.extract.success
    fi
  else
    let BATCH=$BATCH-1
  fi

  echo '#!/bin/bash' > commands.sh
  echo 'if [ ! -e $1.vc.success ];then' >> commands.sh
  echo '  $FREEBAYES -C 2 -0 -O -q 20 -z 0.02 -E 0 -X -u -p 1 -F 0.5 -b $1.alignSorted.bam  -v $1.vcf -f ../$ASM && touch $1.vc.success' >> commands.sh
  echo 'fi' >> commands.sh
  chmod 0755 commands.sh && \
  seq 1 $BATCH |xargs -P $NUM_THREADS -I % ./commands.sh % && \
  for f in $(seq 1 $BATCH);do
    if [ ! -e $f.vc.success ];then
      echo "freebayes failed on batch $f in $BASM.work";
      exit 1
    fi
  done
touch $BASM.vc.success
);
if [ -e ./$BASM.work/$BASM.vc.success ];then
  cat ./$BASM.work/*.vcf > $BASM.vcf
  rm -rf $BASM.work;
  touch $BASM.vc.success
fi
fi

NUMSUB=`grep --text -v '^#'  $BASM.vcf  |perl -ane '{if(length($F[3])==1 && length($F[4])==1){$nerr=1;} print "$F[9]:$nerr\n";}' | awk -F ':' '{if($6>=3 && $4==0) nerr+=$NF}END{print nerr}'`
NUMIND=`grep --text -v '^#' $BASM.vcf  |perl -ane '{if(length($F[3])>1 || length($F[4])>1){$nerr=abs(length($F[3])-length($F[4]));}print "$F[9]:$nerr\n";}' | awk -F ':' '{if($6>=3 && $4==0) nerr+=$NF}END{print nerr}'`
ASMSIZE=`ufasta n50 -S $ASM | awk '{print $2}'`
NUMERR=$(($NUMSUB+$NUMIND))
QUAL=`echo $NUMERR $ASMSIZE | awk '{print 100-$1/$2*100}'`


echo "Substitution Errors: $NUMSUB" > $BASM.report
echo "Insertion/Deletion Errors: $NUMIND" >> $BASM.report
echo "Assembly Size: $ASMSIZE" >> $BASM.report
echo "Consensus Quality: $QUAL" >> $BASM.report
cat $BASM.report

if [ $FIX -gt 0 ];then
  if [ ! -e $BASM.fix.success ];then
    echo "Fixing errors"
    fix_consensus_from_vcf.pl $ASM < $BASM.vcf > $BASM.fixed && touch $BASM.fix.success
  fi
fi

