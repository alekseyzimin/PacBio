#!/bin/bash
refseq=$1;
qryseq=$2;
rnumbases=$3;
nuc_params=$4;
num_cpus=$5


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
echo "NUCmer parameters: $nuc_params"

rm -rf $rundir
mkdir $rundir

split_contig_file.pl $rundir $refseq  $rnumbases 1>/dev/null 2>&1

reffilename=$(basename $refseq)
qryfilename=$(basename $qryseq)

cp $qryseq $rundir/$qryfilename.1

nucdir=$curdir/nucmer.$reffilename.$qryfilename;
tmp_nucdir=$curdir/tmp.nucmer.$reffilename.$qryfilename;

if [ ! -d $nucdir ];then
mkdir -p $tmp_nucdir
(cd $rundir;
parallel -j 1 "if [ ! -e $tmp_nucdir/deltafile.{1}.$qryfilename.1.delta ];then nucmer $nuc_params --delta tmpdeltafile.{1}.$qryfilename.1.delta {1} $qryfilename.1 1>$tmp_nucdir/deltafile.{1}.$qryfilename.1.delta.log 2>&1 && mv tmpdeltafile.{1}.$qryfilename.1.delta $tmp_nucdir/deltafile.{1}.$qryfilename.1.delta;fi" ::: $reffilename.* && mv $tmp_nucdir $nucdir);
fi

if [ -d $nucdir ];then
(cd $rundir;
parallel -j $num_cpus "delta-filter -g -o 20 $nucdir/deltafile.{1}.$qryfilename.1.delta | show-coords -lcHr /dev/stdin > $nucdir/coordsfile.{1}.$qryfilename.coords.tmp && mv $nucdir/coordsfile.{1}.$qryfilename.coords.tmp $nucdir/coordsfile.{1}.$qryfilename.coords" ::: $reffilename.* && cat $nucdir/coordsfile.$reffilename.*.$qryfilename.coords |  awk '{if($4<$5){print $18"/0_"$12" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}else{print $18"/0_"$12" "$19+1" 0 0 0 "$10" "$13-$4+1" "$13-$5+1" "$13" "$1" "$2" "$12" 0"}}' > $curdir/$reffilename.$qryfilename.blasr.out.tmp && mv $curdir/$reffilename.$qryfilename.blasr.out.tmp $curdir/$reffilename.$qryfilename.blasr.out && cd .. && rm -rf $curdir/nucmer.$reffilename.$qryfilename);
fi
