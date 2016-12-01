#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
ESTIMATED_GENOME_SIZE=0
MAX_GAP=2000
MER=15
B=17
d=0.029
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
PB_HC=40;

#parsing arguments
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -m|--masurca_run_path)
    MASURCA_ASSEMBLY_WORK1_PATH="$2"
    shift
    ;;
    -t|--threads)
    NUM_THREADS="$2"
    shift
    ;;
    -M|--alignment_mer)
    MER="$2"
    shift
    ;;
    -B|--alignment_threshold)
    B="$2"
    shift
    ;;
    -D|--density)
    d="$2"
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
KMER=`awk 'BEGIN{min=10000}{if($2<min) min=$2}END{print min}' $KUNITIGLENGTHS`
JF_SIZE=$(stat -c%s $SUPERREADS);
COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}

echo "Running mega-reads correction/assembly"
echo "Using mer size $MER for mapping, B=$B, d=$d"
echo "Using MaSuRCA files from $MASURCA_ASSEMBLY_WORK1_PATH, k-unitig mer $KMER"
echo "Using CA installation from $CA_PATH"
echo "Using $NUM_THREADS threads"
echo "Output prefix $COORDS"

rm -f .rerun
###############removing redundant subreads or reducing the coverage by picking the longest reads##############################
PB_SIZE=$(stat -c%s $PACBIO);
if [ $B -lt 15 ];then
echo "Detected nanopore data";
PACBIO1=$PACBIO;
else
if [ $(($PB_SIZE/$ESTIMATED_GENOME_SIZE)) -gt ${PB_HC} ];then
echo "Pacbio coverage >${PB_HC}x, using ${PB_HC}x of the longest reads";
MAX_GAP=2000
if [ ! -s "pacbio_${PB_HC}xlongest.fa" ] ;then
ufasta extract -f <(ufasta sizes -H $PACBIO | sort -nrk2 -S50% | perl -ane 'BEGIN{$thresh=int("'$ESTIMATED_GENOME_SIZE'")*int("'${PB_HC}'");$n=0}{$n+=$F[1];print $F[0],"\n" if($n<$thresh)}') $PACBIO > pacbio_${PB_HC}xlongest.fa.tmp && mv pacbio_${PB_HC}xlongest.fa.tmp pacbio_${PB_HC}xlongest.fa;
fi
PACBIO1="pacbio_${PB_HC}xlongest.fa";
else
echo "Pacbio coverage <${PB_HC}x, using the longest subreads";
MAX_GAP=1000
if [ ! -s "pacbio_nonredundant.fa" ] ;then
ufasta extract -f <(grep --text '^>' $PACBIO | awk '{print $1}' | awk -F '/' '{split($3,a,"_");print substr($0,2)" "$1"/"$2" "a[2]-a[1]}' | sort -nrk3 -S50% | perl -ane '{if(not(defined($h{$F[1]}))){$h{$F[1]}=1;print $F[0],"\n"}}') $PACBIO > pacbio_nonredundant.fa.tmp && mv pacbio_nonredundant.fa.tmp pacbio_nonredundant.fa;
fi
PACBIO1="pacbio_nonredundant.fa";
fi
fi

if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
echo "Mega-reads pass 1"
if numactl --show 1> /dev/null 2>&1;then
numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $PACBIO1 -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
else
create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $PACBIO1 -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt
fi
touch .rerun
fi

if [ ! -s $COORDS.all_mr.maximal.fa ] || [ -e .rerun ];then
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
perl -ane  '{if($F[0] =~ /^\>/){print substr($F[0],1);}else{ print " ",length($F[0]),"\n";}}' $COORDS.all_mr.fa | sort -nrk2 -S 50%  > $COORDS.mr_sizes.tmp
reduce_sr `wc -l $KUNITIGLENGTHS | perl -ane 'print $F[0]'`  $KUNITIGLENGTHS $KMER $COORDS.mr_sizes.tmp -o $COORDS.reduce.tmp
cat <(awk '{print $1}' $COORDS.reduce.tmp) <(awk '{print $1}'  $COORDS.mr_sizes.tmp) | sort -S 50% |uniq -u > $COORDS.maximal_mr.txt
ufasta extract -f $COORDS.maximal_mr.txt $COORDS.all_mr.fa > $COORDS.all_mr.maximal.fa
rm $COORDS.mr_sizes.tmp $COORDS.reduce.tmp
touch .rerun
fi

if [ ! -s $COORDS.mr.txt ] || [ -e .rerun ];then
echo "Mega-reads pass 2"
if numactl --show 1> /dev/null 2>&1;then
numactl --interleave=all create_mega_reads --stretch-cap 6000 -s $JF_SIZE --psa-min 13 -m $(($MER+2)) -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $(($B-4)) --max-count 2000 -d $d  -r $COORDS.all_mr.maximal.fa  -p $PACBIO1 -o $COORDS.mr.txt.tmp && mv $COORDS.mr.txt.tmp $COORDS.mr.txt
else
create_mega_reads --stretch-cap 6000 -s $JF_SIZE --psa-min 13 -m $(($MER+2)) -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $(($B-4)) --max-count 2000 -d $d  -r $COORDS.all_mr.maximal.fa  -p $PACBIO1 -o $COORDS.mr.txt.tmp && mv $COORDS.mr.txt.tmp $COORDS.mr.txt
fi
touch .rerun
fi

if [ ! -s $COORDS.all.txt ] || [ -e .rerun ];then
echo "Refining alignments"
if which parallel 1> /dev/null 2>&1;then
NUM_PACBIO_READS_PER_BATCH=`grep --text '^>'  $PACBIO1 | wc -l | awk '{bs=int($1/1024);if(bs<1000){bs=1000};if(bs>100000){bs=100000};}END{print bs}'` 
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.mr.txt | add_pb_seq.pl $PACBIO1 | split_matches_file.pl $NUM_PACBIO_READS_PER_BATCH .matches && parallel "refine.sh $COORDS {1} $KMER" ::: .matches.* && cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && rm .matches.* && rm $COORDS.matches*.all.txt.tmp
else
echo "WARNING! GNU parallel not found, proceeding with sequential execution for the refine step"
awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.mr.txt | add_pb_seq.pl $PACBIO1 > .matches.0 && refine.sh $COORDS .matches.0 $KMER && mv $COORDS.matches.0.all.txt.tmp $COORDS.all.txt && rm .matches.0
fi
#awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.mr.txt > $COORDS.all.txt.tmp && mv $COORDS.all.txt.tmp $COORDS.all.txt
touch .rerun
fi

if [ ! -s $COORDS.1.fa ] || [ -e .rerun ];then
echo "Joining"
#allow all joins if high coverage nanopore or pacbio
#if [ $PACBIO1 == $PACBIO  ] || [ $PACBIO1 == "pacbio_${PB_HC}xlongest.fa" ];then
if [ $(($PB_SIZE/$ESTIMATED_GENOME_SIZE)) -gt ${PB_HC} ];then
echo "" > ${COORDS}.1.allowed
else
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
fi
join_mega_reads_trim.onepass.nomatch.pl $PACBIO1 ${COORDS}.1.allowed  $MAX_GAP < ${COORDS}.all.txt 1>$COORDS.1.fa.tmp 2>$COORDS.1.inserts.txt && mv $COORDS.1.fa.tmp $COORDS.1.fa
touch .rerun
fi

if [ ! -s $COORDS.1.frg ] || [ -e .rerun ];then
echo "Generating assembly input files"
make_mr_frg.pl mr 600 < $COORDS.1.fa > $COORDS.1.frg.tmp && mv  $COORDS.1.frg.tmp  $COORDS.1.frg
make_mate_frg.pl < $COORDS.1.fa > $COORDS.1.mates.frg.tmp && mv $COORDS.1.mates.frg.tmp $COORDS.1.mates.frg
rm -rf $CA
fi

TCOVERAGE=20
if [ $ESTIMATED_GENOME_SIZE -gt 1 ];then
MR_SIZE=$(stat -c%s "$COORDS.1.fa");
MCOVERAGE=$((MR_SIZE/ESTIMATED_GENOME_SIZE+1));
if [ $MCOVERAGE -le 5 ];then
echo "Coverage of the mega-reads less than 5 -- using the super reads as well";
SR_FRG=$COORDS.sr.frg
if [ ! -s $SR_FRG ];then
awk '{if($0 ~ /^>/) print $0":super-read"; else print $0}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta | fasta2frg.pl sr 200 > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG;
fi
fi
COVERAGE=`ls $SR_FRG $COORDS.1.frg $COORDS.1.mates.frg $OTHER_FRG 2>/dev/null | xargs stat -c%s | awk '{n+=$1}END{print int(0.85*n/int('$ESTIMATED_GENOME_SIZE'))}'`;
TCOVERAGE=$COVERAGE;
fi

rm -f .rerun
rm -f $CA.log

OVLMIN=`head -n 100000 $SR_FRG $COORDS.1.frg $COORDS.1.mates.frg $OTHER_FRG 2>/dev/null | grep -A 1 '^seq:' |grep -v '^seq:' | grep -v '\-\-' | awk 'BEGIN{minlen=100000}{if(length($1)<minlen && length($1)>=64) minlen=length($1);}END{if(minlen>=250) print "250"; else print minlen-5;}'`
batOptions="-repeatdetect $TCOVERAGE $TCOVERAGE $TCOVERAGE -el $OVLMIN "
OVL_MER=22

echo "Coverage threshold for splitting unitigs is $TCOVERAGE minimum ovl $OVLMIN"

echo "batOptions=$batOptions
cnsConcurrency=$NUM_THREADS
cnsMinFrags=10000
obtMerSize=$OVL_MER
ovlMerSize=$OVL_MER
unitigger=bogart
merylMemory=65536
ovlStoreMemory=65536
utgGraphErrorLimit=1000
utgMergeErrorLimit=1000
utgGraphErrorRate=0.035
utgMergeErrorRate=0.035
ovlCorrBatchSize=100000
ovlCorrConcurrency=4
frgCorrThreads=$NUM_THREADS
mbtThreads=$NUM_THREADS
ovlThreads=2
ovlHashBlockLength=100000000
ovlRefBlockSize=1000000
ovlConcurrency=$NUM_THREADS
doFragmentCorrection=1
doOverlapBasedTrimming=1
doUnitigSplitting=0
doChimeraDetection=normal
merylThreads=$NUM_THREADS
doExtendClearRanges=1
cgwErrorRate=0.15
cgwMergeMissingThreshold=-1
cgwMergeFilterLevel=1
cgwDemoteRBP=0
cnsReuseUnitigs=1" > runCA.spec

echo "Running assembly"
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then 
#need to start from the beginning
runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1 
rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
runCA -s runCA.spec -p genome -d $CA  stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
fi

#at athis point we assume that the unitig consensus is done
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
echo "Unitig consensus failure"
exit;
fi

if [ ! -e "${CA}/deduplicate.success" ]; then
#here we remove overlaps to the reads in duplicate/redundant unitigs and then re-run the unitigger/consensus
deduplicate_unitigs.sh $CA_PATH $CA genome $NUM_THREADS $OVL_MER
runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
runCA -s runCA.spec -p genome -d $CA  stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
echo "Unitig consensus failure after deduplicate"
exit;
else
touch ${CA}/deduplicate.success
fi
fi

#recompute astat if low pacbio coverage
if [ $MCOVERAGE -le 5 ]; then
if [ ! -e ${CA}/recompute_astat.success ];then
recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
fi
fi

#we start from here if the scaffolder has been run or continue here  
runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopAfter=consensusAfterScaffolder $COORDS.1.frg $SR_FRG $COORDS.1.mates.frg $OTHER_FRG 1>> $CA.log 2>&1
rm -rf $CA/8-consensus/*.success $CA/8-consensus/consensus.sh
runCA -s runCA.spec -p genome -d $CA  $COORDS.1.frg $SR_FRG $COORDS.1.mates.frg $OTHER_FRG 1>> $CA.log 2>&1 && echo "Assembly complete. Results are in $CA/9-terminator"

