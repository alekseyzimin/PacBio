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

nucdir=$curdir/nucmer.$reffilename.$qryfilename;
tmp_nucdir=$curdir/tmp.nucmer.$reffilename.$qryfilename;

if [ ! -d $nucdir ];then
mkdir -p $tmp_nucdir
(cd $rundir;
parallel -j $num_cpus "if [ ! -e $tmp_nucdir/deltafile.{1}.{2}.delta ];then nucmer $nuc_params -p tmpdeltafile.{1}.{2} {1} {2} 1>/dev/null 2>&1 && mv tmpdeltafile.{1}.{2}.delta $tmp_nucdir/deltafile.{1}.{2}.delta;fi" ::: $reffilename.* ::: $qryfilename.* && mv $tmp_nucdir $nucdir);
fi

if [ -d $nucdir ];then
(cd $rundir;
for f in $(ls $reffilename.*);do
if [ -e $nucdir/deltafile.$f.$qryfilename.2.delta ];then 
head -n 2 $nucdir/deltafile.$f.$qryfilename.1.delta > $nucdir/deltafile.$f.$qryfilename.merge && tail -n +3 $nucdir/deltafile.$f.$qryfilename.*.delta >> $nucdir/deltafile.$f.$qryfilename.merge
else
ln -s $nucdir/deltafile.$f.$qryfilename.1.delta $nucdir/deltafile.$f.$qryfilename.merge
fi
done
parallel -j $num_cpus "delta-filter -g -o 20 $nucdir/deltafile.{1}.$qryfilename.merge | show-coords -lcHr /dev/stdin > $nucdir/coordsfile.{1}.$qryfilename.coords.tmp && mv $nucdir/coordsfile.{1}.$qryfilename.coords.tmp $nucdir/coordsfile.{1}.$qryfilename.coords" ::: $reffilename.* && cat $nucdir/coordsfile.$reffilename.*.$qryfilename.coords |  awk '{if($4<$5){print $18"/0_"$12" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}else{print $18"/0_"$12" "$19+1" 0 0 0 "$10" "$13-$4+1" "$13-$5+1" "$13" "$1" "$2" "$12" 0"}}' > $curdir/$reffilename.$qryfilename.blasr.out.tmp && mv $curdir/$reffilename.$qryfilename.blasr.out.tmp $curdir/$reffilename.$qryfilename.blasr.out && cd .. && rm -rf $curdir/nucmer.$reffilename.$qryfilename);
fi
