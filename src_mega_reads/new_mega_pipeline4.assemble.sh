#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:~alekseyz/myprogs/masurca-devel/build/inst/bin/:$MYPATH/../build-default/src_jf_aligner:$MYPATH/../src_jf_aligner:$PATH
#arguments
COORDS=$1
KMER=$2
KUNITIGS=$3
SUPERREADS=$4
PACBIO=$5
JF_SIZE=$6

#parameters
MER=15
B=15
d=0.06

#parameters
MER=15
B=15
d=0.06
NUM_THREADS=48

COORDS=$COORDS.$KMER.$MER.$B.$d
if [ -e $COORDS.all.txt ];then
echo "$COORDS.all.txt exists";
exit;
fi

create_mega_reads -s $JF_SIZE -m $MER --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 300 -d $d  -r $SUPERREADS  -p $PACBIO -o $COORDS.txt

perl -ane '{
if($F[0] =~ /^\>/){
$pb=substr($F[0],1);
}else{
$mega_read=$F[8];
@kunis=split(/_/,$mega_read);
$sequence=$F[10];
if(substr($kunis[0],0,-1)>substr($kunis[-1],0,-1)){
        $mega_read=join("_",reverse(@kunis));
        $mega_read=~tr/FR/RF/;
        $sequence=reverse($sequence);
        $sequence=~tr/ACGTNacgtn/TGCANtgcan/;
}
if(not(defined($out{$mega_read}))){
print ">$mega_read\n$sequence\n";
$out{$mega_read}=1;
}
}
}' $COORDS.txt | perl -ane 'BEGIN{$mr_number=0;}{
if($F[0] =~ /^\>/){
$mega_read=substr($F[0],1);
}else{
$sequence=$F[0];
print ">$mr_number\n$sequence\n";
print STDERR "$mega_read\n";
$mr_number++;
@kunis=split(/_/,$mega_read);
$mega_read=join("_",reverse(@kunis));
$mega_read=~tr/FR/RF/;
$sequence=reverse($sequence);
$sequence=~tr/ACGTNacgtn/TGCANtgcan/;
print ">$mr_number\n$sequence\n";
print STDERR "$mega_read\n";
$mr_number++;
}
}' 1>$COORDS.maximal_mr.fa 2>$COORDS.maximal_mr.names

run_big_nucmer_job_parallel.sh pb10x.fasta $COORDS.maximal_mr.fa 1000000 100000000 '--maxmatch -d 0.2 -f -g 200 -l 15 -b 150 -c 100' $NUM_THREADS
show-coords -lcHr  pb10x.fasta.$COORDS.maximal_mr.fa.g.delta | awk '{print $18"/0_"$12" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}' > $COORDS.blasr.out
reconciliate_mega_reads.maximal.nucmer.pl 20 $KMER $COORDS.maximal_mr.fa $COORDS.maximal_mr.names < $COORDS.blasr.out > $COORDS.all.txt

analyze_mega_gaps.sh $COORDS  $KMER > ${COORDS}.allowed; 

eval.sh $COORDS > $COORDS.report &
sleep 20

cd assembly

rm -rf CA.$COORDS.200.[0-9]
rm -rf assembly.CA.$COORDS.200.[0-9]*.*

./commands_iterate.sh ../$COORDS 200 $KMER ../../mega-reads/superreads.frg ../pb10x.fasta CA.$COORDS 2.5

cd ..


