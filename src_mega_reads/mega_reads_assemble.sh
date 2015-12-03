#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
ESTIMATED_GENOME_SIZE=0
#parsing arguments
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -m|--masurca_run_path)
    MASURCA_ASSEMBLY_WORK1_PATH="$2"
    shift
    ;;
    -p|--pacbio)
    PACBIO="$2"
    shift
    ;;
    -a|--assembler_path)
    CA_PATH="$2"
    shift
    ;;
    -e|--estimated-genome-size)
    ESTIMATED_GENOME_SIZE="$2"
    shift
    ;;
    -o|--other_frg)
    OTHER_FRG="$2"
    shift
    ;;
    -v|--verbose)
    set -x
    ;;
    -h|--help|-u|--usage)
    echo "Usage: mega_reads_assemble.sh -m <path to MaSuRCA run work1 folder contents> -p <pacbio reads fasta> -a <path to the location of runCA in wgs-8.2 instalation>"
    exit 0
    ;;
    *)
    echo "Unknown option $1"
    exit 1        # unknown option
    ;;
esac
shift
done

###############checking arguments#########################
if [ ! -e $CA_PATH/runCA ];then
echo "runCA not found at $CA_PATH!";
exit 1;
fi

export PATH=$MYPATH:$CA_PATH:$PATH

if [ ! -e $PACBIO ];then
echo "PacBio reads file $PACBIO not found!";
exit 1 ;
fi

KUNITIGS=$MASURCA_ASSEMBLY_WORK1_PATH/../guillaumeKUnitigsAtLeast32bases_all.fasta
if [ ! -e $KUNITIGS ];then
echo "K-unitigs file $KUNITIGS not found!";
exit 1;
fi

KUNITIGLENGTHS=$MASURCA_ASSEMBLY_WORK1_PATH/kUnitigLengths.txt
if [ ! -e $KUNITIGLENGTHS ];then
echo "K-unitig lengths file $KUNITIGLENGTHS not found!";
exit 1;
fi

SUPERREADS=$MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta
if [ ! -e $SUPERREADS ];then
echo "Super reads file $SUPERREADS not found!";
exit 1;
else
if [ ! -e superReadSequences.named.fasta ];then
perl -ane 'push(@names,$F[0]);
	END{
	open(FILE,"'$SUPERREADS'");
	while($line=<FILE>){
		if($line=~/^>/){
			chomp($line);
			print ">",$names[substr($line,1)],"\n";
		}else{
			print $line;
	}
	}
}' < $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt > superReadSequences.named.fasta.tmp && mv superReadSequences.named.fasta.tmp superReadSequences.named.fasta
fi
if [ -s superReadSequences.named.fasta ];then
SUPERREADS=superReadSequences.named.fasta;
else
echo "Error creating named super-reads file from $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt and $SUPERREADS!";
rm superReadSequences.named.fasta
exit 1;
fi
fi

################setting parameters#########################
MER=15
B=17
d=0.029
KMER=`perl -ane 'BEGIN{$min=10000}{if($F[1]<$min){$min=$F[1]}}END{print $min}' $KUNITIGLENGTHS`
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
REF_BATCH_SIZE=`ls -l $PACBIO | perl -ane '{$s=int($F[4]/10000);$s=100000000 if($s>100000000); print $s;}'`
JF_SIZE=`ls -l $SUPERREADS | perl -ane '{print $F[4]}'`
COORDS=mr.$KMER.$MER.$B.$d
PACBIO_FILE=`basename $PACBIO`;
CA=CA.${COORDS}

echo "Running mega-reads correction/assembly"
echo "Using mer size $MER for mapping, B=$B, d=$d"
echo "Using MaSuRCA files from $MASURCA_ASSEMBLY_WORK1_PATH, k-unitig mer $KMER"
echo "Using CA installation from $CA_PATH"
echo "Using $NUM_THREADS threads"
echo "Output prefix $COORDS"

rm -f .rerun

if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
echo "Mega-reads pass 1"
create_mega_reads -s $JF_SIZE -m $MER --psa-min 14 --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $PACBIO -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
touch .rerun
fi

if [ ! -s $COORDS.all_mr.fa ] || [ -e .rerun ];then
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
}' $COORDS.txt 1> $COORDS.all_mr.fa.tmp && mv $COORDS.all_mr.fa.tmp $COORDS.all_mr.fa 
touch .rerun
fi

if [ ! -s $COORDS.mr.txt ] || [ -e .rerun ];then
echo "Mega-reads pass 2"
create_mega_reads --stretch-cap 6000 -s $JF_SIZE --psa-min 14 -m 17 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B 13 --max-count 1000 -d $d  -r $COORDS.all_mr.fa  -p $PACBIO -o $COORDS.mr.txt.tmp && mv $COORDS.mr.txt.tmp $COORDS.mr.txt
touch .rerun
fi

if [ ! -s $COORDS.all_mr.mr.fa ] || [ -e .rerun ];then
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
}' $COORDS.mr.txt 1> $COORDS.all_mr.mr.fa.tmp  && mv $COORDS.all_mr.mr.fa.tmp $COORDS.all_mr.mr.fa
touch .rerun
fi

if [ ! -s $PACBIO_FILE.$COORDS.maximal_mr.fa.blasr.out ] || [ -e .rerun ];then
echo "Maximal alignment"
perl -ane  '{if($F[0] =~ /^\>/){print substr($F[0],1);}else{ print " ",length($F[0]),"\n";}}' $COORDS.all_mr.mr.fa | sort -nrk2 -S 10%  > $COORDS.mr_sizes.tmp
reduce_sr `wc -l $KUNITIGLENGTHS | perl -ane 'print $F[0]'`  $KUNITIGLENGTHS $KMER $COORDS.mr_sizes.tmp -o $COORDS.reduce.tmp
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
if [ -e .rerun ];then
rm -rf tmp.nucmer.$PACBIO_FILE.$COORDS.maximal_mr.fa;
rm -rf nucmer.$PACBIO_FILE.$COORDS.maximal_mr.fa;
fi
run_big_nucmer_job_parallel_blasr_out.sh $PACBIO $COORDS.maximal_mr.fa $REF_BATCH_SIZE '-d 0.2 -g 200 -l 15 -b 120 -c 100' $NUM_THREADS
touch .rerun
fi

if [ ! -s $COORDS.all.txt ] || [ -e .rerun ];then
echo "Tiling"
reconciliate_mega_reads.maximal.nucmer.pl 20 $KMER $COORDS.maximal_mr.fa $COORDS.maximal_mr.names < $PACBIO_FILE.$COORDS.maximal_mr.fa.blasr.out 1> $COORDS.all.txt.tmp 2>$COORDS.blasr.merged && mv $COORDS.all.txt.tmp $COORDS.all.txt
touch .rerun
fi

if [ ! -s $COORDS.1.fa ] || [ -e .rerun ];then
echo "Joining"
findGapsInCoverageOfPacbios --max-gap-overlap 100  -f $COORDS.blasr.merged > $COORDS.bad_pb.txt.tmp && mv $COORDS.bad_pb.txt.tmp $COORDS.bad_pb.txt
awk 'BEGIN{flag=0}{
        if($0 ~ /^>/){
                flag=0;
                pb=substr($1,2);
        }else{
                flag++;
        }; 
        if(flag>1 && last_mr!=$8){
                l=split(last_mr,a,"_");
                split($8,b,"_");
                k1=int(substr(a[l],1,length(a[l])-1));
                k2=int(substr(b[1],1,length(b[1])-1));
                if(k1<k2){
                        print pb" "$1-$3-last_coord" "k1" "k2
                }else 
                        if(k1>k2){
                                print pb" "$1-$3-last_coord" "k2" "k1
                        }
        };
        last_mr=$8;
        last_coord=$2+$5-$4;
}' ${COORDS}.all.txt |sort -nk3 -k4n -S 20%|uniq -D -f 2 | determineUnjoinablePacbioSubmegas.perl --min-range-proportion 0.15 --min-range-radius 15 > ${COORDS}.1.allowed.tmp && mv ${COORDS}.1.allowed.tmp ${COORDS}.1.allowed
join_mega_reads_trim.onepass.pl $PACBIO ${COORDS}.1.allowed $KMER $COORDS.bad_pb.txt < ${COORDS}.all.txt 1>$COORDS.1.fa.tmp 2>$COORDS.1.inserts.txt && mv $COORDS.1.fa.tmp $COORDS.1.fa
touch .rerun
fi

if [ ! -s $COORDS.1.frg ] || [ ! -s $COORDS.1.mates.frg ] || [ -e .rerun ];then
echo "Generating assembly input files"
make_mr_frg.pl mr 600 < $COORDS.1.fa > $COORDS.1.frg.tmp && mv  $COORDS.1.frg.tmp  $COORDS.1.frg
make_mate_frg.pl < $COORDS.1.fa > $COORDS.1.mates.frg.tmp && mv $COORDS.1.mates.frg.tmp $COORDS.1.mates.frg
rm -rf $CA
fi

if [ $ESTIMATED_GENOME_SIZE -gt 1 ];then
MR_SIZE=$(stat -c%s "$COORDS.1.fa");
COVERAGE=$((MR_SIZE/ESTIMATED_GENOME_SIZE));
if [ $COVERAGE -le 5 ];then
echo "Coverage of the mega-reads less than 5 -- using the super reads as well";
SR_FRG=$COORDS.sr.frg
if [ ! -s $SR_FRG ];then
fasta2frg.pl sr 200 < $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG;
fi
fi
fi
rm -f .rerun

echo "Running assembly"

runCA unitigger=bogart merylMemory=32768 ovlStoreMemory=32768 utgGraphErrorLimit=1000  utgMergeErrorLimit=1000 utgGraphErrorRate=0.04 utgMergeErrorRate=0.04 ovlCorrBatchSize=100000 ovlCorrConcurrency=4 frgCorrThreads=$NUM_THREADS mbtThreads=$NUM_THREADS ovlThreads=2 ovlHashBlockLength=100000000 ovlRefBlockSize=1000000 ovlConcurrency=$NUM_THREADS doFragmentCorrection=1 doOverlapBasedTrimming=1 doUnitigSplitting=0 doChimeraDetection=normal stopAfter=consensusAfterUnitigger consensus=pbutgcns -p genome -d $CA  merylThreads=$NUM_THREADS utgErrorLimit=1000 $COORDS.1.frg $SR_FRG $COORDS.1.mates.frg $OTHER_FRG 1> $CA.log 2>&1

echo "Unitig stats:"
tigStore -g $CA/genome.gkpStore -t $CA/genome.tigStore 2 -U -d sizes 

runCA cnsReuseUnitigs=1 cgwMergeMissingThreshold=-1 cgwMergeFilterLevel=1 cgwDemoteRBP=0 cgwErrorRate=0.25  doFragmentCorrection=1 doOverlapBasedTrimming=1 doUnitigSplitting=0 doChimeraDetection=normal cnsMinFrags=1000 cnsConcurrency=$NUM_THREADS -p genome -d $CA unitigger=bogart merylThreads=$NUM_THREADS utgErrorLimit=1000 $COORDS.1.frg $OTHER_FRG 1>> $CA.log 2>&1

echo "Assembly complete. Results are in $CA/9-terminator"
