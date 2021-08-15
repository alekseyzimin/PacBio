#!/bin/bash
#script to classify reads into parental haplotypes
READSP1=$1
READSP2=$2
MEGAREADS=$3
NUM_THREADS=$4
PID=$$
HAP_THRESH=1.5
set -o pipefail
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PERL5LIB="$MYPATH/../lib/perl/"
export PATH=$MYPATH:$PATH

JF_P1=`basename $READSP1`.jf
JF_P2=`basename $READSP2`.jf
MEGAREADS_P1=`basename $MEGAREADS`.p1
MEGAREADS_P2=`basename $MEGAREADS`.p2
MEGAREADS_BOTH=`basename $MEGAREADS`.both

if [ ! -e count.success ];then
jellyfish count -C -m 31 -t $NUM_THREADS -s 50000000000 -o $JF_P1 -L 3 $READSP1 && \
jellyfish count -C -m 31 -t $NUM_THREADS -s 50000000000 -o $JF_P2 -L 3 $READSP2 && \
touch count.success
fi

if [ ! -e classify.success ];then
cat $JF_P1 $JF_P2 > /dev/null
ufasta split -i $MEGAREADS >($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.1.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.2.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.3.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.4.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.5.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.6.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.7.txt) \
>($MYPATH/classify_reads.pl $JF_P1 $JF_P2 > $PID.8.txt) && \
cat $PID.{1,2,3,4,5,6,7,8}.txt > $PID.counts.txt && \
rm $PID.{1,2,3,4,5,6,7,8}.txt && \
touch classify.success
fi

if [ ! -s extract.success ];then
awk '{if($2>$3){c1=$2+0.0001;c2=$3+0.0001;}else{c2=$2+0.0001;c1=$3+0.0001;}if(c1+c2<20 || c1/c2<='$HAP_THRESH') print $1}' $PID.counts.txt > $MEGAREADS_BOTH.txt && \
awk '{c1=$2+0.0001;c2=$3+0.0001;if(c1+c2>=20 && c1/c2>'$HAP_THRESH') print $1}' $PID.counts.txt > $MEGAREADS_P1.txt && \
awk '{c1=$2+0.0001;c2=$3+0.0001;if(c1+c2>=20 && c2/c1>'$HAP_THRESH') print $1}' $PID.counts.txt > $MEGAREADS_P2.txt && \
ufasta extract -f <(cat $MEGAREADS_BOTH.txt $MEGAREADS_P1.txt ) $MEGAREADS > $MEGAREADS_P1.fa && \
ufasta extract -f <(cat $MEGAREADS_BOTH.txt $MEGAREADS_P2.txt ) $MEGAREADS > $MEGAREADS_P2.fa && \
touch extract.success
fi

#cat p2.fa | /genome2/raid/alekseyz/MaSuRCA/build/inst/bin/make_mr_frg.pl mr 600 > p2.frg 
#cat p1.fa| /genome2/raid/alekseyz/MaSuRCA/build/inst/bin/make_mr_frg.pl mr 600 > p1.frg
