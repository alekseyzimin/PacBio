#!/bin/bash
refseq=$1;
qryseq=$2;
rnumbases=$3;
qnumbases=$4;
nuc_params=$5;
num_cpus=$6


MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

export pid=$$;
export rundir=/dev/shm/$pid.nucmer;
export curdir=$PWD

function cleanup {
if [ $$ -eq $pid ];then
rm -rf $rundir &
fi
}
trap cleanup EXIT SIGHUP SIGINT SIGTERM

echo "Reference: $refseq"
echo "Query: $qryseq"
echo "Reference batch size: $rnumbases"
echo "Query batch size: $qnumbases"
echo "NUCmer parameters: $nuc_params"

rm -rf $rundir
mkdir $rundir

split_contig_file.pl $rundir $refseq  $rnumbases 1>/dev/null 2>&1 &
pid1=$!
split_contig_file.pl $rundir $qryseq  $qnumbases 1>/dev/null 2>&1 &
pid2=$!
wait $pid1 $pid2

reffilename=$(basename $refseq)
qryfilename=$(basename $qryseq)

mkdir -p tmp.nucmer.$reffilename.$qryfilename
(cd $rundir;parallel -j $num_cpus "if [ ! -e $curdir/tmp.nucmer.$reffilename.$qryfilename/deltafile.{1}.{2}.delta ];then nucmer $nuc_params -p tmpdeltafile.{1}.{2} {1} {2} 1>/dev/null 2>&1 && mv tmpdeltafile.{1}.{2}.delta $curdir/tmp.nucmer.$reffilename.$qryfilename/deltafile.{1}.{2}.delta;fi" ::: $reffilename.* ::: $qryfilename.* && mv $curdir/tmp.nucmer.$reffilename.$qryfilename $curdir/nucmer.$reffilename.$qryfilename);

if [ -d nucmer.$reffilename.$qryfilename ];then
cd nucmer.$reffilename.$qryfilename && head -n 2 deltafile.$reffilename.1.$qryfilename.1.delta > $reffilename.$qryfilename.g.delta.tmp && find -type f -name 'deltafile.*.delta'  | xargs cat |grep -v NUCMER | grep -v  nucmer >> $reffilename.$qryfilename.g.delta.tmp && mv $reffilename.$qryfilename.g.delta.tmp ../$reffilename.$qryfilename.g.delta && cd .. && rm -rf nucmer.$reffilename.$qryfilename
fi
