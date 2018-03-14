#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -o pipefail
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
ESTIMATED_GENOME_SIZE=0
MAX_GAP=2000
MER=15
B=17
d=0.029
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
PB_HC=30
KMER=41
#this is the batch size for grid execution
PBATCH_SIZE=2000000000
GRID_ENGINE="SGE"
QUEUE=""
USE_SGE=0
PACBIO=""
NANOPORE=""

function error_exit {
    echo "$1" >&2   
    exit "${2:-1}" 
}


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
        -n|--nanopore)
            NANOPORE="$2"
            shift
            ;;
        -C|--coverage)
            PB_HC="$2"
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
	-G|--use-Grid)
	    USE_SGE="$2"
            shift
	    ;;
        -Pb|--pbatch-size)
            PBATCH_SIZE="$2"
            shift
            ;;
	-q|--queue)
	    QUEUE="$2"
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
echo "Using CABOG from is $CA_PATH"
if [ ! -e $CA_PATH/runCA ];then
    echo "runCA not found at $CA_PATH!";
    exit 1;
fi

export PATH=$MYPATH:$CA_PATH:$PATH

################setting parameters#########################
JF_SIZE=$(stat -L -c%s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta);
if [ $ESTIMATED_GENOME_SIZE -lt 1 ];then 
    echo "Estimated Genome Size is invalid or missing";
    exit;
fi
if [ -s PLOIDY.txt ];then
PLOIDY=`head -n 1 PLOIDY.txt`
else
PLOIDY=$(($JF_SIZE/$ESTIMATED_GENOME_SIZE/2))
fi
if [ $PLOIDY -lt 1 ];then PLOIDY=1; fi
if [ $PLOIDY -gt 2 ];then PLOIDY=2; fi
COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}
echo $CA > CA_dir.txt
echo $PLOIDY > PLOIDY.txt

echo "Running mega-reads correction/assembly"
echo "Using mer size $MER for mapping, B=$B, d=$d"
echo "Estimated Genome Size $ESTIMATED_GENOME_SIZE"
echo "Estimated Ploidy $PLOIDY"
echo "Using $NUM_THREADS threads"
echo "Output prefix $COORDS"

rm -f .rerun
###############removing redundant subreads or reducing the coverage by picking the longest reads##############################

if [ "$PACBIO" = "" ];then
  if [ "$NANOPORE" = "" ];then
    error_exit "must specify either PacBio or Nanopore input reads file with -p or -n"
  else
    LONGREADS=$NANOPORE
  fi
else
  LONGREADS=$PACBIO
fi

if [ ! -e $LONGREADS ];then
  error_exit "Long reads reads file $LONGREADS not found!";
fi

PB_SIZE=$(stat -L -c%s $LONGREADS);
FIRSTCHAR=`zcat -f $LONGREADS | head -c 1`;
if [ "$LONGREADS" = "$NANOPORE" ];then
    if [ ! -s "ont_${PB_HC}xlongest.fa" ] ;then
	echo "Using ${PB_HC}x of the longest ONT reads" 
	if [ "$FIRSTCHAR" = ">" ];then
          zcat -f $LONGREADS |ufasta extract -f <(zcat -f $LONGREADS | ufasta sizes -H | LC_ALL=C sort -nrk2 -S50% | perl -ane 'BEGIN{$thresh=int("'$ESTIMATED_GENOME_SIZE'")*int("'${PB_HC}'")*int("'$PLOIDY'");$n=0}{$n+=$F[1];print $F[0],"\n" if($n<$thresh)}') /dev/stdin > ont_${PB_HC}xlongest.fa.tmp && mv ont_${PB_HC}xlongest.fa.tmp ont_${PB_HC}xlongest.fa || error_exit "failed to extract the best long reads";
        else
	    if [ "$FIRSTCHAR" = "@" ];then
		zcat -f $LONGREADS |fastqToFasta.pl |ufasta extract -f <(zcat -f $LONGREADS | fastqToFasta.pl | ufasta sizes -H | LC_ALL=C sort -nrk2 -S50% | perl -ane 'BEGIN{$thresh=int("'$ESTIMATED_GENOME_SIZE'")*int("'${PB_HC}'")*int("'$PLOIDY'");$n=0}{$n+=$F[1];print $F[0],"\n" if($n<$thresh)}') /dev/stdin > ont_${PB_HC}xlongest.fa.tmp && mv ont_${PB_HC}xlongest.fa.tmp ont_${PB_HC}xlongest.fa || error_exit "failed to extract the best long reads";
	    else
		error_exit "Unknown file format $LONGREADS, exiting";
	    fi
	fi
    fi
    LONGREADS1="ont_${PB_HC}xlongest.fa";
    MAX_GAP=1000
else
    if [ $(($PB_SIZE/$ESTIMATED_GENOME_SIZE/$PLOIDY)) -gt ${PB_HC} ];then
	echo "Pacbio coverage >${PB_HC}x, using ${PB_HC}x of the longest reads";
	MAX_GAP=2000
	if [ ! -s "pacbio_${PB_HC}xlongest.fa" ] ;then
	    if [ "$FIRSTCHAR" = ">" ];then
		zcat -f $LONGREADS |ufasta extract -f <(zcat -f $LONGREADS | grep --text '^>' | awk '{split($1,a,"/");split(a[3],b,"_");len=b[2]-b[1];if($2 ~ /^RQ/){split($2,c,"=");len=int(len*c[2]/0.85);}print substr($1,2)" "len;}'  | LC_ALL=C sort -nrk2 -S50% | perl -ane 'BEGIN{$thresh=int("'$ESTIMATED_GENOME_SIZE'")*int("'${PB_HC}'")*int("'$PLOIDY'");$n=0}{$n+=$F[1];print $F[0],"\n" if($n<$thresh)}') /dev/stdin > pacbio_${PB_HC}xlongest.fa.tmp && mv pacbio_${PB_HC}xlongest.fa.tmp pacbio_${PB_HC}xlongest.fa || error_exit "failed to extract the best long reads";
	    else
		if [ "$FIRSTCHAR" = "@" ];then
		    zcat -f $LONGREADS | fastqToFasta.pl |ufasta extract -f <(zcat -f $LONGREADS | fastqToFasta.pl | grep --text '^>' | awk '{split($1,a,"/");split(a[3],b,"_");len=b[2]-b[1];if($2 ~ /^RQ/){split($2,c,"=");len=int(len*c[2]/0.85);}print substr($1,2)" "len;}'  |LC_ALL=C  sort -nrk2 -S50% | perl -ane 'BEGIN{$thresh=int("'$ESTIMATED_GENOME_SIZE'")*int("'${PB_HC}'")*int("'$PLOIDY'");$n=0}{$n+=$F[1];print $F[0],"\n" if($n<$thresh)}') /dev/stdin  > pacbio_${PB_HC}xlongest.fa.tmp && mv pacbio_${PB_HC}xlongest.fa.tmp pacbio_${PB_HC}xlongest.fa || error_exit "failed to extract the best long reads";
		else
		    error_exit "Unknown file format $LONGREADS, exiting";
		fi
	    fi
	fi
	LONGREADS1="pacbio_${PB_HC}xlongest.fa";
    else
	echo "Pacbio coverage <${PB_HC}x, using the longest subreads";
	MAX_GAP=1000
	if [ ! -s "pacbio_nonredundant.fa" ] ;then
	    if [ "$FIRSTCHAR" = ">" ];then
		zcat -f $LONGREADS |ufasta extract -f <(zcat -f $LONGREADS |grep --text '^>' | awk '{print $1}' | awk -F '/' '{split($3,a,"_");print substr($0,2)" "$1"/"$2" "a[2]-a[1]}' | LC_ALL=C sort -nrk3 -S50% | perl -ane '{if(not(defined($h{$F[1]}))){$h{$F[1]}=1;print $F[0],"\n"}}') /dev/stdin > pacbio_nonredundant.fa.tmp && mv pacbio_nonredundant.fa.tmp pacbio_nonredundant.fa || error_exit "failed to extract the longest subreads"; 
	    else
		if [ "$FIRSTCHAR" = "@" ];then
		    zcat -f $LONGREADS | fastqToFasta.pl |ufasta extract -f <(zcat -f $LONGREADS | fastqToFasta.pl |grep --text '^>' | awk '{print $1}' | awk -F '/' '{split($3,a,"_");print substr($0,2)" "$1"/"$2" "a[2]-a[1]}' | LC_ALL=C sort -nrk3 -S50% | perl -ane '{if(not(defined($h{$F[1]}))){$h{$F[1]}=1;print $F[0],"\n"}}') /dev/stdin > pacbio_nonredundant.fa.tmp && mv pacbio_nonredundant.fa.tmp pacbio_nonredundant.fa || error_exit "failed to extract the longest subreads";
		else
		    error_exit "Unknown file format $LONGREADS, exiting";
		fi
	    fi
	fi
	LONGREADS1="pacbio_nonredundant.fa";
    fi
fi


#first we re-create k-unitigs and super reads with smaller K
if [ ! -e superReadSequences.named.fasta ];then
    echo "Reducing super-read k-mer size"
    awk 'BEGIN{n=0}{if($1~/^>/){}else{print ">sr"n"\n"$0;n+=2;}}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta > superReadSequences.fasta.in.tmp && mv  superReadSequences.fasta.in.tmp  superReadSequences.fasta.in || error_exit "failed to create superReadSequences.fasta.in"
    create_k_unitigs_large_k -q 1 -c $(($KMER-1)) -t $NUM_THREADS -m $KMER -n $(($ESTIMATED_GENOME_SIZE*2)) -l $KMER -f `perl -e 'print 1/'$KMER'/1e5'` $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta.all  | grep --text -v '^>' | perl -ane '{$seq=$F[0]; $F[0]=~tr/ACTGactg/TGACtgac/;$revseq=reverse($F[0]); $h{($seq ge $revseq)?$seq:$revseq}=1;}END{$n=0;foreach $k(keys %h){print ">",$n++," length:",length($k),"\n$k\n"}}' > guillaumeKUnitigsAtLeast32bases_all.fasta.tmp && mv guillaumeKUnitigsAtLeast32bases_all.fasta.tmp guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta || error_exit "failed to create k-unitigs for small k super reads";
    rm -rf work1_mr
    echo "sr 500 50" >>  meanAndStdevByPrefix.pe.txt && createSuperReadsForDirectory.perl -minreadsinsuperread 1 -l $KMER -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta -t $NUM_THREADS -mikedebug work1_mr superReadSequences.fasta.in 1> super1.err 2>&1
    if [ ! -e "work1_mr/superReadSequences.fasta.all" ];then
	error_exit "failed to create super-reads with reduced k-mer size, see super1.err"
    fi
    perl -ane 'push(@names,$F[0]);
END{
  open(FILE,"'work1_mr/superReadSequences.fasta.all'");
  while($line=<FILE>){
    if($line=~/^>/){
      chomp($line);
      print ">",$names[substr($line,1)],"\n";
    }else{
      print $line;
    }
  }
}' < work1_mr/superReadNames.txt > superReadSequences.named.fasta.tmp && mv superReadSequences.named.fasta.tmp superReadSequences.named.fasta || error_exit "failed to create named super-reads file";
    rm superReadSequences.fasta.in
fi

SUPERREADS=superReadSequences.named.fasta
if [ ! -s superReadSequences.named.fasta ];then
    echo "Error creating named super-reads file ";
    exit 1;
fi

KUNITIGS=guillaumeKUnitigsAtLeast32bases_all.$KMER.fasta
if [ ! -s $KUNITIGS ];then
    echo "K-unitigs file $KUNITIGS not found!";
    exit 1;
fi

KUNITIGLENGTHS=work1_mr/kUnitigLengths.txt
if [ ! -s $KUNITIGLENGTHS ];then
    echo "K-unitig lengths file $KUNITIGLENGTHS not found!";
    exit 1;
fi

PBATCHES=$(($(stat -c%s -L $LONGREADS1)/$PBATCH_SIZE));
#if there is one batch then we do not use SGE
if [ $PBATCHES -ge 1001 ];then
    PBATCHES=1000
fi
if [ $PBATCHES -le 1 ];then
    PBATCHES=1
fi

for i in $(seq 1 $PBATCHES);do larr[$i]="lr.batch$i";done;
for i in $(seq 1 $PBATCHES);do lmrOut[$i]="mr.batch$i.txt";done;

if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
    echo "Mega-reads pass 1"

    if [ $PBATCHES -ge 2 ] && [ $USE_SGE -eq 1 ];then
	echo "Running on the grid in $PBATCHES batches";
	if [ "$QUEUE" = "" ];then
	    error_exit "Queue for SGE is undefined, must specify which queue to submit jobs to"
	fi

#here we run on a cluster;first split and then merge alignments

#working inside mr_pass1
	mkdir -p mr_pass1
#running in sub-shell
	( cd mr_pass1;
#split the super-reads
	    if [ ! -e split.success ];then
		ufasta split -i ../$LONGREADS1 ${larr[@]} && touch split.success
	    fi
#creating run scripts
#jf_aligner qsub version
	    echo "#!/bin/sh" > create_mega_reads.sh && \
		echo "if [ ! -e mr.batch\$SGE_TASK_ID.success ];then" >> create_mega_reads.sh && \
		echo "$MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SGE_TASK_ID -o mr.batch\$SGE_TASK_ID.tmp && mv mr.batch\$SGE_TASK_ID.tmp mr.batch\$SGE_TASK_ID.txt && touch mr.batch\$SGE_TASK_ID.success" >> create_mega_reads.sh && \
		echo "else" >> create_mega_reads.sh && \
		echo "echo \"job \$SGE_TASK_ID previously completed successfully\"" >> create_mega_reads.sh && \
		echo "fi"  >> create_mega_reads.sh && chmod 0755 create_mega_reads.sh

#here we use two ways to submit jobs for now. It is straighforward with SGE -- we use the sync option.  For SLURM, we submit the jobs and exit, instructing the used to restart assemble.sh when all jobs finish

#maybe the jobs finished successfully already?
            unset failArr;
            failArr=();
            for i in $(seq 1 $PBATCHES);do
              if [ ! -e mr.batch$i.success ];then
                failArr+=('mr.batch$i')
              fi
            done
            if [ ${#failArr[@]} -ge 1 ];then
              if [ $GRID_ENGINE = "SGE" ];then
              echo "submitting SGE create_mega_reads jobs to the grid"
              qsub -q $QUEUE -cwd -j y -sync y -N "create_mega_reads"  -t 1-$PBATCHES create_mega_reads.sh 1> mqsub2.out 2>&1 || error_exit "create_mega_reads failed on the grid"
              else
              error_exit "submitting SLURM jobs to the grid.  The script will exit now.  Please re-run assemble.sh when all jobs finish."
              fi
            fi

#check if the jobs finished successfully -- needed if we used SGE, redundant for SLURM
            unset failArr;
            failArr=();
            for i in $(seq 1 $PBATCHES);do 
              if [ ! -e mr.batch$i.success ];then
                failArr+=('mr.batch$i')
              fi
            done
            if [ ${#failArr[@]} -ge 1 ];then
              error_exit "${#failArr[@]} create_mega_reads jobs failed in mr_pass1: ${failArr[@]}"
            fi
#cat the results
	    cat ${lmrOut[@]} > ../$COORDS.tmp.txt && mv ../$COORDS.tmp.txt ../$COORDS.txt || error_exit "concatenation of mega-read grid output files failed" ) && rm -rf mr_pass1 || error_exit "mega-reads pass 1 on the grid failed or stopped, please re-run assemble.sh"
    else #single computer
	echo "Running locally in 1 batch";
	if numactl --show 1> /dev/null 2>&1;then
	    numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed";
	else
	    create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed";
	fi
    fi
    touch .rerun
    if  [ ! -s $COORDS.txt ];then
      error_exit "mega-reads pass 1 failed"
    fi
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
}' $COORDS.txt 1> $COORDS.all_mr.fa.tmp && mv $COORDS.all_mr.fa.tmp $COORDS.all_mr.fa || error_exit "failed to extract mega-reads from pass 1 output file";
        ufasta sizes -H $COORDS.all_mr.fa | LC_ALL=C sort -nrk2 -S 50%  > $COORDS.mr_sizes.tmp && \
	reduce_sr `wc -l $KUNITIGLENGTHS | perl -ane 'print $F[0]'`  $KUNITIGLENGTHS $KMER $COORDS.mr_sizes.tmp -o $COORDS.reduce.tmp && \
	cat <(awk '{print $1}' $COORDS.reduce.tmp) <(awk '{print $1}'  $COORDS.mr_sizes.tmp) | LC_ALL=C sort -S 50% | uniq -u > $COORDS.maximal_mr.txt && \
	ufasta extract -f $COORDS.maximal_mr.txt $COORDS.all_mr.fa > $COORDS.all_mr.maximal.fa && \
	rm $COORDS.mr_sizes.tmp $COORDS.reduce.tmp ||  error_exit "failed to create maximal pass 1 mega-reads ";
    touch .rerun
    if  [ ! -s $COORDS.all_mr.maximal.fa ];then
      error_exit "failed to create maximal mega-reads from pass 1"
    fi
fi

#figure out which long reads corrected into one chunk on the first pass
if [ ! -s $COORDS.single.txt ] || [ -e .rerun ];then
    awk 'BEGIN{counter=0}{if($1~ /^>/){if(counter==1){print rn}rn=substr($1,2);counter=0}else{if($8>'$d'*4){counter++}else{counter+=2}}}END{if(counter==1){print rn}}' $COORDS.txt > $COORDS.single.txt.tmp && mv  $COORDS.single.txt.tmp  $COORDS.single.txt || error_exit "failed to extract names of single-chunk mega-reads pass 1";
fi

#here we compute the number of batches to run for secondary create_mega_reads
SBATCHES=$(($(($(($(stat -c%s -L $COORDS.all_mr.maximal.fa)/100000))*$(($(stat -c%s -L $LONGREADS1)/200000))))/$PBATCH_SIZE));

#if fits into 128Gb of RAM, prefer to run on one computer
if [ $(stat -c%s -L $COORDS.all_mr.maximal.fa) -lt 5000000000 ];then
SBATCHES=1
fi

#if there is one batch then we do not use SGE
if [ $SBATCHES -ge 1001 ];then
SBATCHES=1000
fi
if [ $SBATCHES -le 1 ];then
SBATCHES=1
fi
for i in $(seq 1 $SBATCHES);do arr[$i]="sr.batch$i";done;
for i in $(seq 1 $SBATCHES);do arrOut[$i]="coords.batch$i.gz";done;
for i in $(seq 1 $SBATCHES);do mrOut[$i]="mr.batch$i.txt";done;

if [ ! -s $COORDS.mr.txt ] || [ -e .rerun ];then
    echo "Mega-reads pass 2"

    if [ $SBATCHES -ge 2 ] && [ $USE_SGE -eq 1 ];then
      echo "Running on the grid in $SBATCHES batches";
      if [ "$QUEUE" = "" ];then
          error_exit "Queue for SGE is undefined, must specify which queue to submit jobs to"
      fi

#here we run on a cluster;first split and then merge alignments

#working inside mr_pass2
	mkdir -p mr_pass2
#running in sub-shell
	( cd mr_pass2;
#split the super-reads
	    if [ ! -e split.success ];then
		ufasta split -i ../$COORDS.all_mr.maximal.fa ${arr[@]} && touch split.success
	    fi
#creating run scripts
#jf_aligner qsub version
	        echo "#!/bin/bash" > jf_aligner.sh && \
                echo "set -o pipefail" >> jf_aligner.sh && \
		echo "if [ ! -e coords.batch\$SGE_TASK_ID.success ];then" >> jf_aligner.sh && \
		echo "$MYPATH/ufasta extract -v -f ../$COORDS.single.txt ../$LONGREADS1 | $MYPATH/jf_aligner --zero-match -s 1 -m $(($MER+2)) -t $NUM_THREADS -f -B $(($B-4)) --stretch-cap 6000 --max-count $((2000/$SBATCHES)) --psa-min 13 --coords /dev/stdout -u ../$KUNITIGS -k $KMER -H -r sr.batch\$SGE_TASK_ID -p /dev/stdin | ufasta sort -k 2 /dev/stdin | gzip -c -1 > coords.batch\$SGE_TASK_ID.gz && touch coords.batch\$SGE_TASK_ID.success" >> jf_aligner.sh && \
		echo "else" >> jf_aligner.sh && \
		echo "echo \"job \$SGE_TASK_ID previously completed successfully\"" >> jf_aligner.sh && \
		echo "fi"  >> jf_aligner.sh && chmod 0755 jf_aligner.sh

#maybe the jobs finished successfully already?
            unset failArr;
            failArr=();
            for i in $(seq 1 $SBATCHES);do
              if [ ! -e coords.batch$i.success ];then
                failArr+=('coords.batch$i')
              fi
            done
            if [ ${#failArr[@]} -ge 1 ];then
              if [ $GRID_ENGINE = "SGE" ];then
              echo "submitting SGE jf_aligner jobs to the grid"
              qsub -q $QUEUE -cwd -j y -sync y -N "jf_aligner"  -t 1-$SBATCHES jf_aligner.sh 1> jqsub2.out 2>&1 || error_exit "jf_aligner failed on the grid"
              else
              error_exit "submitting SLURM jobs to the grid.  The script will exit now.  Please re-run assemble.sh when all jobs finish."
              fi
            fi
 
#check if the jobs finished successfully -- needed if we used SGE, redundant for SLURM
            unset failArr;
            failArr=();
            for i in $(seq 1 $SBATCHES);do 
              if [ ! -e coords.batch$i.success ];then
                failArr+=('coords.batch$i')
              fi
            done
            if [ ${#failArr[@]} -ge 1 ];then
              error_exit "${#failArr[@]} jf_aligner jobs failed in mr_pass2: ${failArr[@]}"
            fi

#longest path one machine
            echo "#!/bin/bash" > longest_path.sh
            echo "set -o pipefail" >> longest_path.sh
            echo "$MYPATH/merge_coords ${arrOut[@]} |$MYPATH/ufasta extract -v -n \"0\" | ufasta split -i /dev/stdin >($MYPATH/longest_path -t 20  -u ../$KUNITIGS  -k $KMER -d $d -o mr.txt.1.tmp /dev/stdin) >($MYPATH/longest_path -t 20  -u ../$KUNITIGS  -k $KMER -d $d -o mr.txt.2.tmp /dev/stdin) >($MYPATH/longest_path -t 20  -u ../$KUNITIGS  -k $KMER -d $d -o mr.txt.3.tmp /dev/stdin) && cat mr.txt.{1,2,3}.tmp > ../$COORDS.mr.txt" >> longest_path.sh
            chmod 0755 ./longest_path.sh && ./longest_path.sh 
            ) && rm -rf mr_pass2 || error_exit "mega-reads pass 2 on the grid failed or stopped, please re-run assemble.sh"
    else #single computer
        echo "Running locally in 1 batch";
	if numactl --show 1> /dev/null 2>&1;then
	    numactl --interleave=all create_mega_reads --stretch-cap 6000 -s $JF_SIZE --psa-min 13 -m $(($MER+2)) -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $(($B-4)) --max-count 2000 -d $d  -r $COORDS.all_mr.maximal.fa  -p <(ufasta extract -v -f $COORDS.single.txt $LONGREADS1) -o $COORDS.mr.txt.tmp && mv $COORDS.mr.txt.tmp $COORDS.mr.txt || error_exit "mega-reads pass 2 failed";
	else
	    create_mega_reads --stretch-cap 6000 -s $JF_SIZE --psa-min 13 -m $(($MER+2)) -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $(($B-4)) --max-count 2000 -d $d  -r $COORDS.all_mr.maximal.fa  -p <(ufasta extract -v -f $COORDS.single.txt $LONGREADS1) -o $COORDS.mr.txt.tmp && mv $COORDS.mr.txt.tmp $COORDS.mr.txt || error_exit "mega-reads pass 2 failed";
	fi
    fi
    touch .rerun
    if  [ ! -s $COORDS.mr.txt ];then
      error_exit "mega-reads pass 2 failed"
    fi
fi

if [ ! -s $COORDS.all.txt ] || [ -e .rerun ];then
    echo "Refining alignments"
    NUM_LONGREADS_READS_PER_BATCH=`grep --text '^>'  $LONGREADS1 | wc -l | awk '{bs=int($1/1024);if(bs<1000){bs=1000};if(bs>100000){bs=100000};}END{print bs}'` 
    cat <(ufasta extract -f $COORDS.single.txt $COORDS.txt) <(ufasta extract -v -f $COORDS.single.txt $COORDS.mr.txt)| awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' | add_pb_seq.pl $LONGREADS1 | split_matches_file.pl $NUM_LONGREADS_READS_PER_BATCH .matches && ls .matches.* | xargs -P $NUM_THREADS -I % refine.sh $COORDS % $KMER && cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && rm .matches.* && rm $COORDS.matches*.all.txt.tmp 
    touch .rerun
fi

if [ ! -s $COORDS.1.fa ] || [ -e .rerun ];then
    echo "Joining"
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
}' ${COORDS}.all.txt | LC_ALL=C sort -nk3 -k4n -S 20% | determineUnjoinablePacbioSubmegas.perl --min-range-proportion 0.15 --min-range-radius 15 > ${COORDS}.1.allowed.tmp && mv ${COORDS}.1.allowed.tmp ${COORDS}.1.allowed
    join_mega_reads_trim.onepass.nomatch.pl $LONGREADS1 ${COORDS}.1.allowed  $MAX_GAP < ${COORDS}.all.txt 1>$COORDS.1.fa.tmp 2>$COORDS.1.inserts.txt && mv $COORDS.1.fa.tmp $COORDS.1.fa || error_exit "mega-reads joining failed";
    touch .rerun
    if  [ ! -s $COORDS.1.fa ];then
      error_exit "refine/join alignments failed"
    fi
fi

if [ ! -s $COORDS.1.frg ] || [ -e .rerun ];then
    echo "Generating assembly input files"
    awk 'BEGIN{n=0}{if($1~/^>/){}else{print ">sr"n"\n"$0;n+=2;}}'  $COORDS.1.fa  > mr.fa.in && \
	create_k_unitigs_large_k -q 1 -c 30 -t $NUM_THREADS -m 31 -n $ESTIMATED_GENOME_SIZE -l 31 -n $(($ESTIMATED_GENOME_SIZE*2)) -f `perl -e 'print 1/31/1e5'` mr.fa.in   | grep --text -v '^>' | perl -ane '{$seq=$F[0]; $F[0]=~tr/ACTGactg/TGACtgac/;$revseq=reverse($F[0]); $h{($seq ge $revseq)?$seq:$revseq}=1;}END{$n=0;foreach $k(keys %h){print ">",$n++," length:",length($k),"\n$k\n"}}' > guillaumeKUnitigsAtLeast32bases_all.fasta.tmp && mv guillaumeKUnitigsAtLeast32bases_all.fasta.tmp guillaumeKUnitigsAtLeast32bases_all.31.fasta && \
	rm -rf work1_mr1 && \
	createSuperReadsForDirectory.perl -minreadsinsuperread 1 -l 31 -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.31.fasta -t $NUM_THREADS -mikedebug work1_mr1 mr.fa.in  1> super1.err 2>&1 && \
	find_contained_reads.pl work1_mr1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt $COORDS.1.fa > containees.txt && \
	ufasta extract -v -f containees.txt $COORDS.1.fa |make_mr_frg.pl mr 600  > $COORDS.1.frg.tmp && mv  $COORDS.1.frg.tmp  $COORDS.1.frg && \
	make_mate_frg.pl < $COORDS.1.fa > $COORDS.1.mates.frg.tmp && mv $COORDS.1.mates.frg.tmp $COORDS.1.mates.frg && \
        rm -rf $CA work1_mr1 guillaumeKUnitigsAtLeast32bases_all.31.fasta mr.fa.in || error_exit "failed to create mega-reads frg file";
  if  [ ! -s $COORDS.1.frg ];then
    error_exit "failed to create mega-reads frg file"
  fi
fi

TCOVERAGE=20
if [ $ESTIMATED_GENOME_SIZE -gt 1 ];then
    MR_SIZE=$(stat -c%s -L "$COORDS.1.fa");
    MCOVERAGE=$(($MR_SIZE/$ESTIMATED_GENOME_SIZE/$PLOIDY+1));
    if [ $MCOVERAGE -le 5 ];then
	echo "Coverage of the mega-reads less than 5 -- using the super reads as well";
	SR_FRG=$COORDS.sr.frg
	if [ ! -s $SR_FRG ];then
	    awk '{if($0 ~ /^>/) print $0":super-read"; else print $0}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta | fasta2frg.pl sr 200 > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG || error_exit "failed to create super-reads frg file";
	fi
    fi
    COVERAGE=`ls $SR_FRG $COORDS.1.frg $COORDS.1.mates.frg $OTHER_FRG 2>/dev/null | xargs stat -c%s | awk '{n+=$1}END{cov=int(n/int('$ESTIMATED_GENOME_SIZE')/int('$PLOIDY')); if(cov<15) cov=15; print cov;}'`;
    TCOVERAGE=$COVERAGE;
fi

rm -f .rerun
rm -f $CA.log

OVLMIN=`head -n 100000 $SR_FRG $COORDS.1.frg $COORDS.1.mates.frg $OTHER_FRG 2>/dev/null | grep -A 1 '^seq:' |grep -v '^seq:' | grep -v '\-\-' | awk 'BEGIN{minlen=100000}{if(length($1)<minlen && length($1)>=64) minlen=length($1);}END{if(minlen>=250) print "250"; else print minlen-1;}'`
batOptions="-repeatdetect $TCOVERAGE $TCOVERAGE $TCOVERAGE -el $OVLMIN "
OVL_MER=22

echo "Coverage threshold for splitting unitigs is $TCOVERAGE minimum ovl $OVLMIN"
let NUM_THREADSd4=$(($NUM_THREADS/4+1))
if [ $USE_SGE -ge 1 ];then
OVL_THREADS=4
else
OVL_THREADS=2
fi

echo "batOptions=$batOptions
useGrid=$USE_SGE
cnsOnGrid=0
gridSubmitCommand=qsub
cnsConcurrency=$NUM_THREADS
cnsMinFrags=10000
obtMerSize=$OVL_MER
ovlMerSize=$OVL_MER
unitigger=bogart
merylMemory=65536
ovlStoreMemory=65536
utgGraphErrorLimit=1000
utgMergeErrorLimit=1000
utgGraphErrorRate=0.03
utgMergeErrorRate=0.03
ovlCorrBatchSize=100000
ovlCorrConcurrency=6
frgCorrThreads=$NUM_THREADS
frgCorrConcurrency=6
mbtThreads=$NUM_THREADS
ovlThreads=$OVL_THREADS
ovlHashBlockLength=10000000
ovlRefBlockSize=1000000
ovlConcurrency=$NUM_THREADS
doOverlapBasedTrimming=1
doUnitigSplitting=0
doChimeraDetection=normal
merylThreads=$NUM_THREADS
stoneLevel=0
doExtendClearRanges=0
computeInsertSize=0
maxRepeatLength=12000
cgwErrorRate=0.1
cgwMergeMissingThreshold=-1
cgwMergeFilterLevel=1
cgwDemoteRBP=0
cnsReuseUnitigs=1" > runCA.spec


echo "Running assembly"
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then 
  #need to start from the beginning
  #this is helpful for re-starts
  rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh
    $CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopAfter=consensusAfterUnitigger $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
    #this is a fix for sometimes failing fragment correction
    if [ ! -e "${CA}/4-unitigger/unitigger.success" ]; then
      rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh
      $CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopAfter=consensusAfterUnitigger $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG doFragmentCorrection=0 1>> $CA.log 2>&1
    fi
    rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
    $CA_PATH/runCA -s runCA.spec -p genome -d $CA  stopAfter=consensusAfterUnitigger $COORDS.1.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
fi

#at athis point we check if the unitig consensus is done
if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
    echo "Assembly stopped or failed, see $CA.log" 
    exit;
fi

if [ ! -e "${CA}/deduplicate.success" ]; then
#here we remove overlaps to the reads in duplicate/redundant unitigs and then re-run the unitigger/consensus
    deduplicate_unitigs.sh $CA_PATH $CA genome $NUM_THREADS $OVL_MER $PLOIDY
fi

if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
  #after deduplicate we need to rebuild the unitigs, we rerun CA on deduplicated overlapStore
    $CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopAfter=consensusAfterUnitigger $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
    rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
    $CA_PATH/runCA -s runCA.spec -p genome -d $CA  stopAfter=consensusAfterUnitigger cnsConcurrency=$(($NUM_THREADS/2+1)) $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
fi

if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
    echo "Assembly stopped or failed, see $CA.log"
    exit;
fi

#recompute astat if low pacbio coverage
if [ $MCOVERAGE -le 5 ]; then
    if [ ! -e ${CA}/recompute_astat.success ];then
	echo "Recomputing A-stat"
	recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
	touch ${CA}/recompute_astat.success
    fi
fi

#we start from here if the scaffolder has been run or continue here  
$CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopAfter=consensusAfterScaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
rm -rf $CA/8-consensus/*.success $CA/8-consensus/consensus.sh
$CA_PATH/runCA -s runCA.spec -p genome -d $CA  cnsConcurrency=$(($NUM_THREADS/2+1)) $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1 && \
    echo "Mega-reads initial assembly complete." || echo "Assembly stopped or failed, see $CA.log"

