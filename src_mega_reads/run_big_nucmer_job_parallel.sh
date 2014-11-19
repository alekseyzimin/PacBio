#!/bin/bash
refseq=$1;
qryseq=$2;
rnumbases=$3;
qnumbases=$4;
nuc_params=$5;
num_cpus=$6

export pid=$$;
export rundir=/dev/shm/$refseq.$qryseq.$pid.nucmer;

function cleanup {
if [ $$ -eq $pid ];then
rm -rf $rundir;
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

(cd $rundir;parallel -j $num_cpus "nucmer $nuc_params -p deltafile.{1}.{2} {1} {2} 1>/dev/null 2>&1" ::: $refseq.* ::: $qryseq.* ;)

head -n 2 $rundir/deltafile.$refseq.1.$qryseq.1.delta > $refseq.$qryseq.g.delta
cat  $rundir/deltafile*.delta |grep -v NUCMER | grep -v  $refseq >> $refseq.$qryseq.g.delta


