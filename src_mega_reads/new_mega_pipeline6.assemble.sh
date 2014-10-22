#!/bin/bash
#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:~alekseyz/myprogs/masurca-devel/build/inst/bin/:$MYPATH/../build-default:$MYPATH/../build-default/src_jf_aligner:$MYPATH/../build-default/src_mega_reads:$PATH
#arguments
COORDS=$1
KMER=$2
KUNITIGS=$3
KUNITIGLENGTHS=$4
SUPERREADS=$5
PACBIO=$6
JF_SIZE=$7

#parameters
MER=13
B=25
d=0.1
NUM_THREADS=16

COORDS=$COORDS.$KMER.$MER.$B.$d
if [ ! -e $COORDS.all.txt ];then

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

create_mega_reads -s $JF_SIZE -m $MER --max-match -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 300 -d $d  -r $COORDS.all_mr.fa  -p $PACBIO -o $COORDS.mr.txt

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
reduce_sr 148481 $KUNITIGLENGTHS $KMER $COORDS.mr_sizes.tmp -o $COORDS.reduce.tmp
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

fi

show-coords -lcHr  pb10x.fasta.$COORDS.maximal_mr.fa.g.delta | perl -e '{@lines=();$pb="";while($line=<STDIN>){chomp($line);$line=~s/^\s+//;@f=split(/\s+/,$line);if(not($pb eq $f[17])){if($#lines>=0){print join("\n",sort by18thenby0 @lines),"\n";}$pb=$f[17];@lines=();}push(@lines,$line);} sub by18thenby0 { my @ffa=split(/\s+/,$a); my @ffb=split(/\s+/,$b);return($ffa[18] <=> $ffb[18] || $ffa[0] <=> $ffb[0])}}' | merge_matches_coords_file.pl 2000 | awk '{print $18" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}' > $COORDS.blasr.merged
findGapsInCoverageOfPacbios --max-gap-overlap `perl -e '{print int('$KMER'*1.2)}'`  -f $COORDS.blasr.merged > $COORDS.bad_pb.txt

analyze_mega_gaps.sh $COORDS $KMER | determineUnjoinablePacbioSubmegas.perl --min-range-proportion 0.15 --min-range-radius 15 > ${COORDS}.allowed

#eval.sh $COORDS $KMER > $COORDS.report & 

cd assembly

rm -rf CA.$COORDS.200.[0-9]*
rm -rf assembly.CA.$COORDS.200.[0-9]*

./commands_iterate.sh ../$COORDS $KMER ../superreads.frg ../$PACBIO CA.$COORDS 5

cd ..


