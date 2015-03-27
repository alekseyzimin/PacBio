#!/bin/bash
#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:/home/alekseyz/myprogs/masurca-devel/build/inst/bin/:$MYPATH/../build-default:$MYPATH/../build-default/src_jf_aligner:$MYPATH/../build-default/src_mega_reads:/home/alekseyz/wgs-8.2/Linux-amd64/bin/:$PATH
#arguments
COORDS=$1
KMER=$2
KUNITIGS=$3
KUNITIGLENGTHS=$4
SUPERREADS=$5
PACBIO=$6
JF_SIZE=$7

#parameters
MER=15
B=17
d=0.04
NUM_THREADS=48

COORDS=$COORDS.$KMER.$MER.$B.$d
CA=CA.${COORDS}

if [ ! -e $COORDS.blasr.out ];then

create_mega_reads -s $JF_SIZE -m $MER --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 300 -d $d  -r $SUPERREADS  -p $PACBIO -o $COORDS.txt

#first we reduce to maximal
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
}' $COORDS.txt 1> $COORDS.all_mr.fa 

create_mega_reads -s $JF_SIZE -m $MER --max-match -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 300 -d 0.05  -r $COORDS.all_mr.fa  -p $PACBIO -o $COORDS.mr.txt

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
}' $COORDS.mr.txt 1> $COORDS.all_mr.mr.fa 


perl -ane  '{if($F[0] =~ /^\>/){print substr($F[0],1);}else{ print " ",length($F[0]),"\n";}}' $COORDS.all_mr.mr.fa | sort -nrk2 -S 10%  > $COORDS.mr_sizes.tmp
reduce_sr `wc -l $4 | perl -ane 'print $F[0]'`  $KUNITIGLENGTHS $KMER $COORDS.mr_sizes.tmp -o $COORDS.reduce.tmp
cat <(awk '{print $1}' $COORDS.reduce.tmp) <(awk '{print $1}'  $COORDS.mr_sizes.tmp) | sort -S 10% |uniq -u > $COORDS.maximal_mr.txt
extractreads.pl $COORDS.maximal_mr.txt $COORDS.all_mr.mr.fa 1 | perl -ane 'BEGIN{$mr_number=0;}{
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
print STDERR "$mega_read\n";
$mr_number++;
}
}' 1>$COORDS.maximal_mr.fa 2>$COORDS.maximal_mr.names

run_big_nucmer_job_parallel.sh $PACBIO $COORDS.maximal_mr.fa 2000000 200000000 '-d 0.2 -g 200 -l 15 -b 120 -c 100' $NUM_THREADS
delta-filter -g -o 20 $PACBIO.$COORDS.maximal_mr.fa.g.delta > $PACBIO.$COORDS.maximal_mr.fa.gg.delta
show-coords -lcHr  $PACBIO.$COORDS.maximal_mr.fa.gg.delta | awk '{if($4<$5){print $18"/0_"$12" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}else{print $18"/0_"$12" "$19+1" 0 0 0 "$10" "$13-$4+1" "$13-$5+1" "$13" "$1" "$2" "$12" 0"}}' > $COORDS.blasr.out

reconciliate_mega_reads.maximal.nucmer.pl 20 $KMER $COORDS.maximal_mr.fa $COORDS.maximal_mr.names < $COORDS.blasr.out 1> $COORDS.all.txt 2>$COORDS.blasr.merged

findGapsInCoverageOfPacbios --max-gap-overlap 100  -f $COORDS.blasr.merged > $COORDS.bad_pb.txt

analyze_mega_gaps.sh $COORDS  $KMER | determineUnjoinablePacbioSubmegas.perl --min-range-proportion 0.15 --min-range-radius 15 > ${COORDS}.1.allowed

fi

join_mega_reads_trim.onepass.pl $PACBIO ${COORDS}.1.allowed $KMER $COORDS.bad_pb.txt < ${COORDS}.all.txt > $COORDS.1.fa;

fasta2frg.pl mr 600 < $COORDS.1.fa > $COORDS.1.frg;
#perl -ane 'BEGIN{$n=0}{if($F[0]=~/^>/){print ">$n\n";$n++}else{print "$F[0]\n"}}' $SUPERREADS | fasta2frg.pl sr  > $COORDS.sr.frg;

make_mate_frg.pl < $COORDS.1.fa > $COORDS.1.mates.frg

runCA unitigger=bogart  merylMemory=8192 utgGraphErrorLimit=1000  utgMergeErrorLimit=1000 utgGraphErrorRate=0.04 utgMergeErrorRate=0.04 ovlCorrBatchSize=5000 ovlCorrConcurrency=$NUM_THREADS frgCorrThreads=$NUM_THREADS mbtThreads=$NUM_THREADS ovlThreads=2 ovlHashBlockLength=100000000 ovlRefBlockSize=1000 ovlConcurrency=$NUM_THREADS doFragmentCorrection=1 doOverlapBasedTrimming=1 doUnitigSplitting=0 doChimeraDetection=normal doUnitigSplitting=0 stopAfter=unitigger cnsMinFrags=500 cnsConcurrency=32 -p genome -d $CA  merylThreads=$NUM_THREADS utgErrorLimit=1000 $COORDS.1.frg    $COORDS.1.mates.frg  1> $CA.log 2>&1

echo "preliminary stats:"
tigStore -g $CA/genome.gkpStore -t $CA/genome.tigStore 2 -U -d sizes -s 12000000

runCA cnsReuseUnitigs=1 cgwMergeMissingThreshold=-1 cgwMergeFilterLevel=1 cgwDemoteRBP=0 cgwErrorRate=0.25 merylMemory=8192 utgGraphErrorLimit=1000 utgGraphErrorRate=0.035 utgMergeErrorLimit=1000 utgMergeErrorRate=0.045 ovlCorrBatchSize=5000 ovlCorrConcurrency=$NUM_THREADS frgCorrThreads=$NUM_THREADS mbtThreads=$NUM_THREADS ovlThreads=2 ovlHashBlockLength=100000000 ovlRefBlockSize=10000 ovlConcurrency=$NUM_THREADS doFragmentCorrection=1 doOverlapBasedTrimming=1 doUnitigSplitting=0 doChimeraDetection=normal doUnitigSplitting=0 cnsMinFrags=300 cnsConcurrency=32 -p genome -d $CA unitigger=bogart merylThreads=$NUM_THREADS utgErrorLimit=1000 $COORDS.1.frg  1>> $CA.log 2>&1

#split_long_unitigs.pl cr < ${CA}/9-terminator/genome.ctg.fasta 2>assembly.${CA}.short_contigs.fa | fasta2frg.pl cr > assembly.${CA}.contig_reads.frg

#runCA merylMemory=8192 unitigger=bogart utgGraphErrorLimit=1000 utgGraphErrorRate=0.0 utgMergeErrorLimit=1000 utgMergeErrorRate=0.01 ovlMerThreshold=5 ovlMinLen=1000 doFragmentCorrection=0 doUnitigSplitting=0 doChimeraDetection=off stopAfter=unitigger cnsMinFrags=100 cnsConcurrency=16 -p genome -d ${CA}u ovlThreads=$NUM_THREADS merylThreads=$NUM_THREADS doOverlapBasedTrimming=0 utgErrorLimit=100000 assembly.${CA}.contig_reads.frg  1>> $CA.log 2>&1

#tigStore -g ${CA}u/genome.gkpStore -t ${CA}u/genome.tigStore 2 -U -d sizes -s 12000000

#final output
#tigStore -g ${CA}u/genome.gkpStore -t ${CA}u/genome.tigStore  3 -U -d consensus >assembly.${CA}u.fa

