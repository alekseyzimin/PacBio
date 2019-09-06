#!/bin/bash

set -o pipefail
set -e
NUM_THREADS=4
MEM=16000000000

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            export NUM_THREADS="$2"
            shift
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
            echo "Usage:  evaluate_consensus_error_rate.sh -a <assembly contigs or scaffolds> -r <Illumina reads fastq> -t <number of threads>"
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

if [ ! -e $BASM.index.success ];then 
rm -f $BASM.map.success
$BWA index $ASM -p $BASM.bwa && touch $BASM.index.success
fi

if [ ! -e $BASM.map.success ];then
rm -f $BASM.sort.success
$BWA mem -SP -t 64 $BASM.bwa $READS 2>bwasterr |samtools view -bhS /dev/stdin > $BASM.unSorted.bam && touch $BASM.map.success
fi
#here we are doing sorting, indexing and variant calling in parallel, per input contig/scaffold
if [ ! -e $BASM.call.success ];then
grep '^>' $BASM |awk '{print substr($1,2)}' > $BASM.names
mkdir -p $BASM.work
(cd $BASM.work
for f in $(cat ../$BASM.names);do
  $SAMTOOLS view -h ../$BASM.unSorted.bam $f |$SAMTOOLS view -S -b /dev/stdin > $f.unSorted.bam
done
echo '#!/bin/bash' > commands.sh
echo 'if [ ! -e $1.sort.success]'
echo '  $SAMTOOLS sort $1.unSorted.bam $1.alignSorted && $SAMTOOLS index $1.alignSorted.bam && touch  $1.sort.success'
echo 'fi'
echo 'if [ ! -e $1.vc.success]'
echo '  $FREEBAYES -C 2 -0 -O -q 20 -z 0.02 -E 0 -X -u -p 1 -F 0.5 -b $1.alignSorted.bam  -v $1.vcf -f ../$ASM && touch $1.vc.success'
echo 'fi'
chmod 0755 commands.sh
cat ../$BASM.names |xargs -P $NUM_THREADS -I % ./commands.sh
cat *.vcf > ../$BASM.vcf);
touch $BASM.call.success
fi



#if [ ! -e $BASM.sort.success ];then
#rm -f $BASM.vc.success
#$SAMTOOLS sort -m $MEM  $BASM.unSorted.bam $BASM.alignSorted && \
#$SAMTOOLS index $BASM.alignSorted.bam && touch  $BASM.sort.success
#fi

#if [ ! -e $BASM.vc.success ];then
#rm -f $BASM.calc.success
#freebayes -C 2 -0 -O -q 20 -z 0.02 -E 0 -X -u -p 1 -F 0.5 -b $BASM.alignSorted.bam  -v $BASM.vcf -f $ASM && touch $BASM.vc.success
#fi

NUMSUB=`grep -v "#" $BASM.vcf  |perl -ane '{if(length($F[3])==1 && length($F[4])==1){$nerr=1;} print "$F[9]:$nerr\n";}' | awk -F ':' '{if($6>=3 && $4==0) nerr+=$NF}END{print nerr}'`
NUMIND=`grep -v "#" $BASM.vcf  |perl -ane '{if(length($F[3])>1 || length($F[4])>1){$nerr=abs(length($F[3])-length($F[4]));}print "$F[9]:$nerr\n";}' | awk -F ':' '{if($6>=3 && $4==0) nerr+=$NF}END{print nerr}'`
ASMSIZE=`ufasta n50 -S $ASM | awk '{print $2}'`
NUMERR=$(($NUMSUB+$NUMIND))
QUAL=`echo $NUMERR $ASMSIZE | awk '{print 100-$1/$2*100}'`


echo "Substitution Errors: $NUMSUB" > $BASM.report
echo "Insertion/Deletion Errors: $NUMIND" >> $BASM.report
echo "Assembly Size: $ASMSIZE" >> $BASM.report
echo "Consensus Quality: $QUAL" >> $BASM.report

