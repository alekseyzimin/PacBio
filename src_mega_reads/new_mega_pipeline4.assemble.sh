#!/bin/bash
COORDS=$1
KMER=$2
MER=$3
B=$4
d=$5


COORDS=$COORDS.$KMER.$MER.$B.$d
if [ -e $COORDS.all.txt ];then
echo "$COORDS.all.txt exists";
exit;
fi

~/myprogs/PacBio/src_jf_aligner/create_mega_reads -s 20000000 -m $MER -k 70 -u ../assembly_k70/guillaumeKUnitigsAtLeast32bases_all.fasta -t 48 -B $B --max-count 300 -d $d  -r ../assembly_k70/work1/superReadSequences.named.fasta  -p pb10x.fasta -o $COORDS.txt

perl -ane 'BEGIN{$mrn=0;}{
if($F[0] =~ /^\>/){
$pb=substr($F[0],1);
}else{
$mega_read=$F[8];
$sequence=$F[10];
if(not(defined($out{$mega_read}))){
print ">$mrn\n$sequence\n";
$out{$mega_read}=$mrn;
$mrn++;
}
$pacbios{"$mega_read:$out{$mega_read}"}.="$pb:$F[6] ";
}
}END{
foreach $m(keys %pacbios){
print STDERR "$m $pacbios{$m}\n";
}}' $COORDS.txt 1> $COORDS.all_mr.mr.fa 2>$COORDS.mr.pb.txt

./run_big_nucmer_job_parallel.sh pb10x.fasta $COORDS.all_mr.mr.fa 1000000 100000000 '--maxmatch -d 0.2 -f -g 200 -l 15 -b 150' 48
show-coords -lcHr  pb10x.fasta.$COORDS.all_mr.mr.fa.g.delta | awk '{print $18"/0_"$12" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}' > $COORDS.blasr.out
~/myprogs/PacBio/src_mega_reads/reconciliate_mega_reads.nucmer.pl 20 $KMER $COORDS.all_mr.mr.fa $COORDS.mr.pb.txt< $COORDS.blasr.out > $COORDS.all.txt

./analyze_mega_gaps.sh $COORDS  $KMER > ${COORDS}.allowed; 

./eval.sh $COORDS > $COORDS.report &
sleep 20

cd assembly

rm -rf CA.$COORDS.200.[0-9]
rm -rf assembly.CA.$COORDS.200.[0-9]*.*

./commands_iterate.sh ../$COORDS 200 $KMER ../../mega-reads/superreads.frg ../pb10x.fasta CA.$COORDS 2.5

cd ..


