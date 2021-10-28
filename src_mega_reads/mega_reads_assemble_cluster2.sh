#!/bin/bash
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/bin/bash
set -o pipefail
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export MYPATH
ESTIMATED_GENOME_SIZE=0
#minimum 50
MAX_GAP=1000
MER=17
B=15
d=0.02
NUM_THREADS=`cat /proc/cpuinfo |grep ^processor |wc -l`
PB_HC=30
KMER=0
#this is the batch size for grid execution
PBATCH_SIZE=5000000000
GRID_ENGINE="SGE"
QUEUE="all.q"
USE_GRID=0
PACBIO=""
NANOPORE=""
NANOPORE_RNA=0
ONEPASS=0
OVLMIN_DEFAULT=250
FLYE_ASSEMBLY=0
MASURCA_ASSEMBLY_WORK1_PATH="work1"
POSTFIX=""
MIN_PROPORTION="0.25"
MIN_RADIUS="400"
SUPERREADS=superReadSequences.named.fasta
KUNITIGS=guillaumeKUnitigsAtLeast32bases_all.fasta
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

trap abort 1 2 15
function abort {
log "Aborted"
kill 0
exit 1
}


function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2   
    exit "${2:-1}" 
}

USAGE="Usage: 
mega_reads_assemble_cluster2.sh <arguments>
arguments:

MANDATORY
-k int      MANDATORY k-mer size used to construct k-unitigs
-a string   MANDATORY path to the contig assembler bin folder
-e int      MANDATORY estimated genome size

-p string   input PacBio reads 
or 
-n string   input Nanopore reads

OPTIONAL, have defaults
-m string   path to MaSuRCA run work1 folder default:work1
-u string   k-unitigs file name default:guillaumeKUnitigsAtLeast32bases_all.fasta
-t int      number of threads default:number of CPUs
-M int      k-mer size for alignments to long reads default:17
-B double   alignment threshold % bases default:15.0
-D double   density of k-mers default:0.02
-C int      max coverage of PacBio/nanopore reads to use default:25
-G bool     whether to use the grid 0 is no | 1 is yes default:0
-E string   grid engine SGE or SLURM default:SGE
-Pb int     batch size in number of bases for grid execution default:5000000000
-q string   grid queue to use default:all.q
-F          use Flye assembler for building contigs default:not set
-o string   other FRG file to include in assembly, ignored if -F is set
-v          verbose default:not set
-h|-u       this message
" 

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
        -k|--kunitig_mer)
            KMER="$2"
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
        -r|--rnaseq)
            NANOPORE="$2"
            NANOPORE_RNA=1
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
	-G|--use-grid)
	    USE_GRID="$2"
            shift
	    ;;
        -E|--grid-engine)
            GRID_ENGINE="$2"
            shift
            ;;
        -u|--kunitigs)
            KUNTIGS="$2"
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
        -F|--Flye)
            FLYE_ASSEMBLY=1;
            ;;
	-v|--verbose)
	    set -x
	    ;;    
	-h|--help|-u|--usage)
	    echo "$USAGE"
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
if [ $FLYE_ASSEMBLY -gt 0 ];then
    log "Using Flye from $CA_PATH"
    if [ ! -e $CA_PATH/flye ];then
	error_exit "flye not found at $CA_PATH!";
    fi
else
    log "Using CABOG from $CA_PATH"
    if [ ! -e $CA_PATH/runCA ];then
	error_exit "runCA not found at $CA_PATH!";
    fi
fi
export PATH=$MYPATH:$CA_PATH:$PATH

if [ ! -s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta ];then
  error_exit "super reads file not found or size zero, you can try deleting $MASURCA_ASSEMBLY_WORK1_PATH folder and re-generating assemble.sh, also check if guillaumeKUnitigsAtLeast32bases_all.fasta is not empty";
fi

if [ ! -s $KUNITIGS ];then
  error_exit "k-unitigs file $KUNITIGS not found or size zero (-u)"
fi

if [ $KMER -lt 4 ];then
  error_exit "improper value for k-unitig k-mer (-k)"
fi

if [ $MER -lt 4 ];then
  error_exit "improper value for alignment mer (-M)"
fi

################setting parameters#########################
JF_SIZE=$(stat -L -c%s $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta);
if [ $ESTIMATED_GENOME_SIZE -lt 1 ];then 
    error_exit "Estimated Genome Size is invalid or missing";
fi
if [ -s PLOIDY.txt ];then
    PLOIDY=`head -n 1 PLOIDY.txt| awk '{if($1<=1) print "1"; else if($1>=2) print "2"; else print "1"}'`
else
    PLOIDY=$(($JF_SIZE/$ESTIMATED_GENOME_SIZE/2))
fi
if [ $PLOIDY -lt 1 ];then PLOIDY=1; fi
if [ $PLOIDY -gt 2 ];then PLOIDY=2; fi
echo $PLOIDY > PLOIDY.txt

#we need to increase MER for very large genomes 
if [ $ESTIMATED_GENOME_SIZE -gt 10000000000 ];then
  MER=18
  B=12
fi



COORDS=mr.$KMER.$MER.$B.$d
FLYE=flye.${COORDS}
CA=CA.${COORDS}
echo $CA > CA_dir.txt
if [ $FLYE_ASSEMBLY -lt 1 ];then
echo $CA > CA_dir.txt
else
echo $FLYE > FLYE_dir.txt
fi


log "Running mega-reads correction/assembly"
log "Using mer size $MER for mapping, B=$B, d=$d"
log "Estimated Genome Size $ESTIMATED_GENOME_SIZE"
log "Estimated Ploidy $PLOIDY"
log "Using $NUM_THREADS threads"
log "Output prefix $COORDS"

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

LR_SIZE=$(stat -L -c%s $LONGREADS);
FIRSTCHAR=`zcat -f $LONGREADS | head -c 1`
LR_COVERAGE=$(($LR_SIZE/$ESTIMATED_GENOME_SIZE/$PLOIDY))

if [ $LR_COVERAGE -gt $PB_HC ];then
    OVLMIN_DEFAULT=499
fi

LONGREADS1=longest_reads.${PB_HC}x.fa
if [ ! -s $LONGREADS1 ];then
#here we do the initial pre-correction of the long reads and pick the best ones to use for the remaining steps
  $MYPATH/correct_with_k_unitigs.sh -k 19  -t $NUM_THREADS -e $ESTIMATED_GENOME_SIZE -i pe.cor.fa -l $LONGREADS -c $PB_HC -o $LONGREADS1
  if [ ! -s $LONGREADS1 ];then
    error_exit "Failed to pre-correct $LONGREADS file, please check your data!"
  fi
fi

#first we re-create k-unitigs and super reads with smaller K
if [ ! -s $SUPERREADS ];then
    perl -ane 'push(@names,$F[0]);
END{
  open(FILE,"'$MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta'");
  while($line=<FILE>){
    if($line=~/^>/){
      chomp($line);
      print ">",$names[substr($line,1)],"\n";
    }else{
      print $line;
    }
  }
}' < $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt > $SUPERREADS.tmp && mv $SUPERREADS.tmp $SUPERREADS || error_exit "failed to create named super-reads file";
fi

if [ ! -s $SUPERREADS ];then
    error_exit "Error creating named super-reads file ";
fi

if [ ! -s $KUNITIGS ];then
    error_exit "K-unitigs file $KUNITIGS not found; check your inputs"
fi

KUNITIGLENGTHS=$MASURCA_ASSEMBLY_WORK1_PATH/kUnitigLengths.txt
if [ ! -s $KUNITIGLENGTHS ];then
    error_exit "K-unitig lengths file $KUNITIGLENGTHS not found!";
fi

PBATCHES=$(($(stat -c%s -L $LONGREADS1)/$PBATCH_SIZE));
#if there is one batch then we do not use SGE
if [ $PBATCHES -ge 512 ];then
    PBATCHES=512
fi
if [ $PBATCHES -le 1 ];then
    PBATCHES=1
fi

for i in $(seq 1 $PBATCHES);do larr[$i]="lr.batch$i";done;
for i in $(seq 1 $PBATCHES);do lmrOut[$i]="mr.batch$i.txt";done;

if [ ! -s $COORDS.txt ] || [ -e .rerun ];then
    log "Computing mega-reads"

    if [ $PBATCHES -ge 2 ] && [ $USE_GRID -eq 1 ];then
	log "Running on the grid in $PBATCHES batches";
	if [ "$QUEUE" = "" ];then
	    error_exit "Queue for the grid is undefined, must specify which queue to submit jobs to"
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
#creating run scripts for SGE or SLURM
            if [ $GRID_ENGINE = "SGE" ];then
		echo "#!/bin/sh" > create_mega_reads.sh && \
		    echo "if [ ! -e mr.batch\$SGE_TASK_ID.success ];then" >> create_mega_reads.sh && \
                    echo "if numactl --show 1> /dev/null 2>&1;then" >> create_mega_reads.sh && \
		    echo "numactl --interleave=all $MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 11  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SGE_TASK_ID -o mr.batch\$SGE_TASK_ID.tmp && mv mr.batch\$SGE_TASK_ID.tmp mr.batch\$SGE_TASK_ID.txt && touch mr.batch\$SGE_TASK_ID.success" >> create_mega_reads.sh && \
                    echo "else" >> create_mega_reads.sh && \
                    echo "$MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 11  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SGE_TASK_ID -o mr.batch\$SGE_TASK_ID.tmp && mv mr.batch\$SGE_TASK_ID.tmp mr.batch\$SGE_TASK_ID.txt && touch mr.batch\$SGE_TASK_ID.success" >> create_mega_reads.sh && \
                    echo "fi"  >> create_mega_reads.sh
		    echo "else" >> create_mega_reads.sh && \
		    echo "echo \"job \$SGE_TASK_ID previously completed successfully\"" >> create_mega_reads.sh && \
		    echo "fi"  >> create_mega_reads.sh && chmod 0755 create_mega_reads.sh
            else
		echo "#!/bin/sh" > create_mega_reads.sh && \
		    echo "if [ ! -e mr.batch\$SLURM_ARRAY_TASK_ID.success ];then" >> create_mega_reads.sh && \
                    echo "if numactl --show 1> /dev/null 2>&1;then" >> create_mega_reads.sh && \
                    echo "numactl --interleave=all $MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 11  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SLURM_ARRAY_TASK_ID -o mr.batch\$SLURM_ARRAY_TASK_ID.tmp && mv mr.batch\$SLURM_ARRAY_TASK_ID.tmp mr.batch\$SLURM_ARRAY_TASK_ID.txt && touch mr.batch\$SLURM_ARRAY_TASK_ID.success" >> create_mega_reads.sh && \
                    echo "else" >> create_mega_reads.sh && \
		    echo "$MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 11  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SLURM_ARRAY_TASK_ID -o mr.batch\$SLURM_ARRAY_TASK_ID.tmp && mv mr.batch\$SLURM_ARRAY_TASK_ID.tmp mr.batch\$SLURM_ARRAY_TASK_ID.txt && touch mr.batch\$SLURM_ARRAY_TASK_ID.success" >> create_mega_reads.sh && \
                    echo "fi"  >> create_mega_reads.sh
		    echo "else" >> create_mega_reads.sh && \
		    echo "echo \"job \$SLURM_ARRAY_TASK_ID previously completed successfully\"" >> create_mega_reads.sh && \
		    echo "fi"  >> create_mega_reads.sh && chmod 0755 create_mega_reads.sh
            fi

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
		    log "submitting SGE create_mega_reads jobs to the grid"
		    qsub -q $QUEUE -cwd -j y -sync y -N "create_mega_reads"  -t 1-$PBATCHES create_mega_reads.sh 1> mqsub2.out 2>&1
                    #we submit the jobs the second time to re-run (hopefully successfully) any jobs that may hav failed
                    qsub -q $QUEUE -cwd -j y -sync y -N "create_mega_reads"  -t 1-$PBATCHES create_mega_reads.sh 1> mqsub2.out 2>&1 || error_exit "create_mega_reads failed on the grid"
		else
		    echo " "
		    echo "To submit SLURM jobs, please run"
		    echo " "
		    echo "sbatch -D `pwd` -J create_mega_reads -a 1-$PBATCHES -n $NUM_THREADS -p $QUEUE -N 1 mr_pass1/create_mega_reads.sh"
		    echo " "
		    echo "Please re-run assemble.sh when all jobs finish. If you get this message again, it means that some jobs failed, simply re-submit again using the above command."
		    echo " "
		    exit 1
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
		error_exit "${#failArr[@]} create_mega_reads jobs failed in mr_pass1: ${failArr[@]}, re-run assemble.sh"
            fi
#cat the results
	    cat ${lmrOut[@]} > ../$COORDS.tmp.txt && mv ../$COORDS.tmp.txt ../$COORDS.txt || error_exit "concatenation of mega-read grid output files failed" ) && rm -rf mr_pass1 || error_exit "mega-reads pass 1 on the grid exited, please re-run assemble.sh"
    else #single computer
	log "Running locally in 1 batch";
        #if previous temporary file exists -- run failure, then we continue after deleting the last entry in the $COORDS.txt.tmp; cannot use ufasta because the file may be corrupted/incomplete
        if [ -e $COORDS.txt.tmp ];then #if found previous temp file
          if [ ! -s $COORDS.txt.tmp ];then # and it is empty
            if [ -s $COORDS.txt.done ];then #but we already tried to continue
              mv $COORDS.txt.done $COORDS.txt.tmp #then move the previous done file to the temp file to not lose the work
            fi
          fi
        fi
        if [ -s $COORDS.txt.tmp ];then # found failed run, need to be careful below, the file may be corrupted
        log "Found $COORDS.txt.tmp file, attempting to continue from where the previous run left off"
          cat -n $COORDS.txt.tmp | grep '>' > $COORDS.txt.headers.tmp && \
          mv $COORDS.txt.headers.tmp $COORDS.txt.headers && \
          if [ -s $COORDS.txt.headers ];then
            head -n `tail -n 1 $COORDS.txt.headers | perl -ane '{print $F[0]-1}'` $COORDS.txt.tmp > $COORDS.txt.done.tmp && \
            mv $COORDS.txt.done.tmp $COORDS.txt.done
          else
            echo "1 >_\n2 >_" > $COORDS.txt.headers
            touch $COORDS.txt.done
          fi
          if numactl --show 1> /dev/null 2>&1;then #if numactl then interleave memory
            numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p <(ufasta extract -v -f <(head -n -1 $COORDS.txt.headers | awk '{print substr($2,2)}') $LONGREADS1)  -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt 
          else
            create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p <(ufasta extract -v -f <(head -n -1 $COORDS.txt.headers | awk '{print substr($2,2)}') $LONGREADS1)  -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt
          fi
          rm -f $COORDS.txt.headers
          if [ -e $COORDS.txt ];then # if success of mega-reads
            cat $COORDS.txt.done $COORDS.txt > $COORDS.txt.tmp1 && rm $COORDS.txt.done && mv $COORDS.txt.tmp1 $COORDS.txt || error_exit "mega-reads pass 1 failed, likely ran out of disk space" 
          else #mega-reads failed again
            cat $COORDS.txt.done $COORDS.txt.tmp > $COORDS.txt.tmp1 && mv $COORDS.txt.tmp1 $COORDS.txt.tmp && rm $COORDS.txt.done
            error_exit "mega-reads pass 1 failed, please re-generate and re-run assemble.sh to continue"
          fi
        else #no failed run
          if numactl --show 1> /dev/null 2>&1;then #if numactl then interleave memory
            numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed"
          else
            create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed"
          fi
        fi #failed run
    fi #single computer
    touch .rerun
    if  [ ! -s $COORDS.txt ];then
	error_exit "mega-reads pass 1 failed"
    fi
fi

if [ ! -s ${COORDS}.all.txt ] || [ -e .rerun ];then
    log "Refining alignments"
    NUM_LONGREADS_READS_PER_BATCH=`grep --text '^>'  $LONGREADS1 | wc -l | awk '{bs=int($1/1024);if(bs<1000){bs=1000};if(bs>100000){bs=100000};}END{print bs}'` 
    awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' $COORDS.txt | \
    $MYPATH/add_pb_seq.pl $LONGREADS1 | \
    $MYPATH/split_matches_file.pl $NUM_LONGREADS_READS_PER_BATCH .matches && \
    ls .matches.* | xargs -P $NUM_THREADS -I % $MYPATH/refine.sh $COORDS % $KMER && \
    cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && \
    rm .matches.* && rm $COORDS.matches*.all.txt.tmp 
    touch .rerun
fi

if [ $NANOPORE_RNA -gt 0 ];then
    join_mega_reads_trim.onepass.ref.pl <$COORDS.all.txt 1>$COORDS.transcripts.fa 2>$COORDS.transcripts.err
    exit
fi

if [ ! -s ${COORDS}.1$POSTFIX.allowed ] || [ -e .rerun ];then
    log "Computing allowed merges"
    cat ${COORDS}.all.txt | $MYPATH/determineUnjoinablePacbioSubmegas.perl --min-range-proportion $MIN_PROPORTION --min-range-radius $MIN_RADIUS > ${COORDS}.1$POSTFIX.allowed.tmp && mv ${COORDS}.1$POSTFIX.allowed.tmp ${COORDS}.1$POSTFIX.allowed || error_exit "computing allowed merges failed"
    touch .rerun
fi

if [ ! -e ${COORDS}.1$POSTFIX.unjoined.fa ] || [ ! -e ${COORDS}.1$POSTFIX.to_join.fa ] || [ -e .rerun ];then
    log "Joining"
    ALL_SIZE=$(stat -L -c%s ${COORDS}.all.txt);
    if [ $ALL_SIZE -gt 10000000 ];then
#use 4 processes
	ufasta split -i <($MYPATH/add_pb_seq.pl $LONGREADS1 < ${COORDS}.all.txt) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.1 2>$COORDS.1$POSTFIX.to_join.fa.tmp.1) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.2 2>$COORDS.1$POSTFIX.to_join.fa.tmp.2) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.3 2>$COORDS.1$POSTFIX.to_join.fa.tmp.3) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.4 2>$COORDS.1$POSTFIX.to_join.fa.tmp.4) && \
            cat $COORDS.1$POSTFIX.fa.tmp.{1,2,3,4} | awk '{if($1!~/^>/ && $1~/>/){split($1,a,">");print a[1];if(a[2]!="") print ">"a[2];}else{print $0}}' > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.unjoined.fa && \
            cat $COORDS.1$POSTFIX.to_join.fa.tmp.{1,2,3,4} | awk '{if($1!~/^>/ && $1~/>/){split($1,a,">");print a[1];if(a[2]!="") print ">"a[2];}else{print $0}}' > $COORDS.1$POSTFIX.to_join.fa.tmp && mv $COORDS.1$POSTFIX.to_join.fa.tmp $COORDS.1$POSTFIX.to_join.fa && \
            rm $COORDS.1$POSTFIX.fa.tmp.{1,2,3,4} $COORDS.1$POSTFIX.to_join.fa.tmp.{1,2,3,4} || error_exit "mega-reads joining failed"
    else
#use 1 process
	$MYPATH/add_pb_seq.pl $LONGREADS1 < ${COORDS}.all.txt | $MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp 2>$COORDS.1$POSTFIX.to_join.fa.tmp && mv $COORDS.1$POSTFIX.to_join.fa.tmp $COORDS.1$POSTFIX.to_join.fa && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.unjoined.fa || error_exit "mega-reads joining failed"
    fi
    touch .rerun
fi

if [ ! -s $COORDS.1$POSTFIX.fa ] || [ -e .rerun ];then
    log "Gap consensus"    
#making consensus for large gaps
    mkdir -p ${COORDS}.join_consensus.tmp && \

#start subshell execution

    (cd ${COORDS}.join_consensus.tmp;
	TOJOIN_BATCHES=$(($(stat -c%s -L ../${COORDS}.1$POSTFIX.to_join.fa)/100000000))
	if [ $TOJOIN_BATCHES -le 1 ];then
	    TOJOIN_BATCHES=1
	fi
	if [ $TOJOIN_BATCHES -ge 500 ];then
	    TOJOIN_BATCHES=500
	fi
	for i in $(seq 1 $TOJOIN_BATCHES);do ref_names[$i]="ref.$i.fa";done;
	for i in $(seq 1 $TOJOIN_BATCHES);do merges_names[$i]="merges.$i.txt";done;

	awk 'BEGIN{flag=1}{if($2>int("'$MAX_GAP'")*.75 && $6==1) {if($3==prev3 && $4==prev4) flag++; else flag=1;  print $5-$2" "$1" "$3" "$4;prev3=$3;prev4=$4}}'  ../${COORDS}.1$POSTFIX.allowed  > qrys.txt && \
	    perl -ane '{$index="$F[2] $F[3]";if(not(defined($hl{$index})) || $ho{$index}>$F[0]){$hl{$index}=join(" ",@F);$ho{$index}=$F[0];}}END{foreach $k(keys %hl){print $hl{$k},"\n";}}' qrys.txt > refs.txt && \
	    if [ ! -s qrys.all.fa ]; then $MYPATH/ufasta extract -f <(awk '{print $2}' qrys.txt) ../$LONGREADS1 > qrys.all.fa.tmp && mv qrys.all.fa.tmp qrys.all.fa; fi && \
	    $MYPATH/ufasta extract -v -f <(awk '{print $2}' refs.txt) qrys.all.fa | $MYPATH/ufasta one > qrys.fa && \
	    $MYPATH/ufasta extract -f <(awk '{print $2}' refs.txt) qrys.all.fa | $MYPATH/ufasta one > refs.fa && \
	    if [ -s qrys.fa ] && [ -s refs.fa ];then #if no joins to make
	    perl -ane '{
       $h{$F[1]}="$F[2]_$F[3]";
       }END{
       $flag=0;
       open(FILE,"refs.fa");
       while($line=<FILE>){
         if($line=~/^>/){
           chomp($line);
           @f=split(/\s+/,$line);
           if(defined($h{substr($f[0],1)})){
             $line=<FILE>;
             chomp($line);
             $hseq{$h{substr($f[0],1)}}=$line;
           }
         }
       }
       foreach $name(keys %hseq){
         print ">$name\n$hseq{$name}\n";
       }}' refs.txt > refs.renamed.fa && \
           rm -f ${ref_names[@]} && $MYPATH/ufasta split -i  refs.renamed.fa ${ref_names[@]} && \
           $MYPATH/split_reads_to_join.pl qrys.txt to_blasr ${ref_names[@]} < qrys.fa && \
           perl -ane '{if($F[0] =~ /^>/){$rn=$F[0];}else{$seq=$F[0]; $seq=~ tr/a-zA-Z//s; print "$rn\n$F[0]\n" if(length($seq)>length($F[0])*0.1);}}' ../${COORDS}.1$POSTFIX.to_join.fa | $MYPATH/split_reads_to_join.pl qrys.txt to_join ${ref_names[@]} && \
           grep '^>' --text ../$COORDS.1$POSTFIX.to_join.fa | perl -ane '{($rn,$coord)=split(/\./,substr($F[0],1));$h{$rn}.=substr($F[0],1)." ";}END{foreach $r(keys %h){@f=split(/\s+/,$h{$r}); for ($i=0;$i<$#f;$i++){print $f[$i]," ",$f[$i+1],"\n"}}}' > valid_join_pairs.txt && \
           for F in $(seq 1 $TOJOIN_BATCHES);do echo ">_0" >> to_join.$F.fa;echo "ACGT" >> to_join.$F.fa;done && \
           for F in $(seq 1 $TOJOIN_BATCHES);do echo ">_0" >> to_blasr.$F.fa;echo "ACGT" >> to_blasr.$F.fa;done 


#only try blasr/pbdagcon for one pass
           if [ ! -e pass_one ];then
            touch pass_one
#gap consensus on the grid, if set
#try to do gap consensus with blasr; some jobs may fail
            echo "#!/bin/bash" > ./do_consensus.sh
		echo "set -o pipefail" >> ./do_consensus.sh
                if [ $USE_GRID -eq 1 ]; then
                  if [ $GRID_ENGINE = "SGE" ];then
                    echo "TASK_ID=\$SGE_TASK_ID"  >> ./do_consensus.sh
                    log "Using SGE grid"
                  else
                    echo "TASK_ID=\$SLURM_ARRAY_TASK_ID" >> ./do_consensus.sh
                    log "Using SLURM grid queue $QUEUE"
                  fi
                else
                  echo "TASK_ID=\$1" >> ./do_consensus.sh
                fi
		echo "if [ ! -e consensus.\$TASK_ID.success ];then" >> ./do_consensus.sh
                if [ $USE_GRID -eq 1 ]; then
		  echo "$MYPATH/../CA8/Linux-amd64/bin/blasr to_blasr.\$TASK_ID.fa   ref.\$TASK_ID.fa  -minMatch 15 -nproc $NUM_THREADS -bestn 10 -m 5 2>blasr.err | \\" >> ./do_consensus.sh && \
                  echo "sort -k6 -S5% | $MYPATH/../CA8/Linux-amd64/bin/pbdagcon -j $NUM_THREADS -t 0 -c 1 /dev/stdin  2>pbdagcon.err | awk -F 'N' '{if(\$1 == \"\") print \"ACGT\"; else print \$1}' > join_consensus.\$TASK_ID.fasta && \\" >> ./do_consensus.sh
                  echo "$MYPATH/nucmer --delta /dev/stdout --batch 10000 -l 17 -c 51 -L 200 -t $NUM_THREADS to_join.\$TASK_ID.fa join_consensus.\$TASK_ID.fasta 2>/dev/null | \\" >> ./do_consensus.sh
                else
                  echo "$MYPATH/../CA8/Linux-amd64/bin/blasr to_blasr.\$TASK_ID.fa   ref.\$TASK_ID.fa  -minMatch 15 -nproc 16 -bestn 10 -m 5 2>blasr.err | \\" >> ./do_consensus.sh 
                  echo "sort -k6 -S5% | $MYPATH/../CA8/Linux-amd64/bin/pbdagcon -j 8 -t 0 -c 1 /dev/stdin  2>pbdagcon.err | awk -F 'N' '{if(\$1 == \"\") print \"ACGT\"; else print \$1}' > join_consensus.\$TASK_ID.fasta && \\" >> ./do_consensus.sh
                  echo "$MYPATH/nucmer --delta /dev/stdout --batch 10000 -l 17 -c 51 -L 200 -t 16 to_join.\$TASK_ID.fa join_consensus.\$TASK_ID.fasta 2>/dev/null | \\" >> ./do_consensus.sh
                fi
                echo "$MYPATH/filter_delta_file_for_qrys.pl qrys.txt | \\" >> ./do_consensus.sh && \
                echo "$MYPATH/show-coords -lcHq -I 88 /dev/stdin > coords.\$TASK_ID && cat coords.\$TASK_ID | \\" >> ./do_consensus.sh && \
                echo "$MYPATH/extract_merges_mega-reads.pl join_consensus.\$TASK_ID.fasta  valid_join_pairs.txt > merges.\$TASK_ID.txt && touch consensus.\$TASK_ID.success" >> ./do_consensus.sh && \
                echo "fi" >> ./do_consensus.sh && \
                echo "exit 0" >> ./do_consensus.sh && \
		chmod 0755 ./do_consensus.sh

                if [ $USE_GRID -eq 1 ]; then
                  if [ $GRID_ENGINE = "SGE" ];then
                    qsub -cwd -j y -sync y -N "join_mega_reads"  -t 1-$TOJOIN_BATCHES do_consensus.sh 1> cqsub2.out 2>&1 || error_exit "join consensus failed on the grid"  
                  else
                        echo " "
                        echo "To submit SLURM jobs, please run"
                        echo " "
                        echo "(cd ${COORDS}.join_consensus.tmp && sbatch -D `pwd` -J join_mega_reads -a 1-$TOJOIN_BATCHES -n $NUM_THREADS -p $QUEUE -N 1 do_consensus.sh);"
                        echo " "
                        echo "Please re-run assemble.sh when all jobs finish. If you get this message again, it means that some jobs failed, simply re-submit again using the above command."
                        echo " "
                        exit 0 
                  fi
                else
                  seq 1 $TOJOIN_BATCHES | xargs -P 4 -I % ./do_consensus.sh %
                fi #use_grid

            fi #pass one

#re-do the failed jobs with flye -- more reliable but less sensisitive
            echo "#!/bin/bash" > ./do_consensus.sh && \
                echo "set -o pipefail" >> ./do_consensus.sh && \
                echo "if [ ! -e consensus.\$1.success ];then" >> ./do_consensus.sh && \
                echo "$MYPATH/../Flye/bin/flye --polish-target ref.\$1.fa --iterations 1 --pacbio-raw to_blasr.\$1.fa --out-dir \$1.polish.tmp --threads \$2 2>flye.\$1.err && \\" >> ./do_consensus.sh && \
                echo "ufasta one \$1.polish.tmp/polished_1.fasta | \\" >> ./do_consensus.sh && \
                echo "awk '{if(\$1 ~ /^>/) header=\$1; else print header\"/0_\"length(\$1)\"\\n\"\$1}' | tee join_consensus.\$1.fasta | \\" >> ./do_consensus.sh && \
                echo "$MYPATH/nucmer --delta /dev/stdout --batch 10000 -l 17 -c 51 -L 200 -t \$2 to_join.\$1.fa /dev/stdin 2>/dev/null | \\" >> ./do_consensus.sh && \
                echo "$MYPATH/filter_delta_file_for_qrys.pl qrys.txt | \\" >> ./do_consensus.sh && \
                echo "$MYPATH/show-coords -lcHq -I 88 /dev/stdin > coords.\$1 && cat coords.\$1 | \\" >> ./do_consensus.sh && \
                echo "$MYPATH/extract_merges_mega-reads.pl join_consensus.\$1.fasta  valid_join_pairs.txt > merges.\$1.txt && touch consensus.\$1.success" >> ./do_consensus.sh && \
                echo "fi" >> ./do_consensus.sh && \
                echo "exit 0" >> ./do_consensus.sh && \
                chmod 0755 ./do_consensus.sh
 

            if [ $TOJOIN_BATCHES -le 1 ];then #use all CPUs for one batch
              ./do_consensus.sh 1 $NUM_THREADS
            else #do 4 batches at a time for efficiency
              seq 1 $TOJOIN_BATCHES | xargs -P 4 -I % ./do_consensus.sh % $(($NUM_THREADS/4+1))
            fi

            cat merges.[0-9]*.txt |perl -ane '{if($F[2] eq "F"){$merge="$F[0] $F[3]";}else{$merge="$F[3] $F[0]";} if(not(defined($h{$merge}))|| $h{$merge} > $F[1]+$F[4]){$hl{$merge}=join(" ",@F);$h{$merge}=$F[1]+$F[4];}}END{foreach $k(keys %hl){print $hl{$k},"\n"}}' > merges.best.txt && \
		$MYPATH/merge_mega-reads.pl < merges.best.txt | \
		$MYPATH/create_merged_mega-reads.pl ../${COORDS}.1$POSTFIX.to_join.fa merges.best.txt > ${COORDS}.1$POSTFIX.joined.fa.tmp && \
		mv ${COORDS}.1$POSTFIX.joined.fa.tmp  ../${COORDS}.1$POSTFIX.joined.fa && \
		cat ${merges_names[@]} > /dev/null && \
		touch join_consensus.success
                else #if no joins to make
                cp ../${COORDS}.1$POSTFIX.to_join.fa ../${COORDS}.1$POSTFIX.joined.fa.tmp && mv ../${COORDS}.1$POSTFIX.joined.fa.tmp ../${COORDS}.1$POSTFIX.joined.fa && touch join_consensus.success
                fi)

#end subshell execution

    if [ -e ${COORDS}.join_consensus.tmp/join_consensus.success ];then
        cat ${COORDS}.1$POSTFIX.joined.fa ${COORDS}.1$POSTFIX.unjoined.fa  > ${COORDS}.1$POSTFIX.fa.tmp && \
            mv ${COORDS}.1$POSTFIX.fa.tmp ${COORDS}.1$POSTFIX.fa && rm -rf ${COORDS}.join_consensus.tmp 
    else
        log "Warning! Some or all gap consensus jobs failed, see files in ${COORDS}.join_consensus.tmp, however this is fine and assembly can proceed normally"
        if [ -s ${COORDS}.1$POSTFIX.joined.fa ];then
          if [ $FLYE_ASSEMBLY -gt 0 ];then
            cat $COORDS.1$POSTFIX.joined.fa $COORDS.1$POSTFIX.unjoined.fa | $MYPATH/ufasta extract -f <($MYPATH/ufasta sizes -H <(cat $COORDS.1$POSTFIX.joined.fa $COORDS.1$POSTFIX.unjoined.fa) | sort -nrk2,2 -S 10% | awk '{n+=$2;if(n<549755813888) print $1;}') > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.fa
          else
            cat $COORDS.1$POSTFIX.joined.fa $COORDS.1$POSTFIX.unjoined.fa  > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.fa
          fi
        else
          if [ $FLYE_ASSEMBLY -gt 0 ];then
            cat $COORDS.1$POSTFIX.joined.fa $COORDS.1$POSTFIX.to_join.fa | $MYPATH/ufasta extract -f <($MYPATH/ufasta sizes -H <(cat $COORDS.1$POSTFIX.to_join.fa $COORDS.1$POSTFIX.unjoined.fa) | sort -nrk2,2 -S 10% | awk '{n+=$2;if(n<549755813888) print $1;}') > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.fa
          else
            cat $COORDS.1$POSTFIX.unjoined.fa $COORDS.1$POSTFIX.to_join.fa  > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.fa
          fi
        fi
    fi
    touch .rerun
    if  [ ! -s $COORDS.1$POSTFIX.fa ];then
        error_exit "Gap consensus failed"
    fi
fi

if [ $FLYE_ASSEMBLY -gt 0 ];then
    if [ ! -s "$FLYE/assembly.fasta" ];then
        log "Running assembly with Flye"
	$CA_PATH/flye -t $NUM_THREADS --nano-corr $COORDS.1$POSTFIX.fa -g $ESTIMATED_GENOME_SIZE --kmer-size 21 -m 2500 -o $FLYE -i 0 1>flye.log 2>&1
    fi
    (cd $FLYE && $MYPATH/samba.sh -r assembly.fasta -q $LONGREADS -t $NUM_THREADS && mv assembly.fasta.scaffolds.fa assembly.scaffolds.fasta);
else
    if [ ! -s $COORDS.1.frg ] || [ ! -s $COORDS.1.mates.frg ] || [ -e .rerun ];then
	log "Generating assembly input files"
	awk 'BEGIN{n=0}{if($1~/^>/){}else{print ">sr"n"\n"$0;n+=2;}}'  $COORDS.1$POSTFIX.fa  > mr.fa.in && \
	    create_k_unitigs_large_k -q 1 -c 30 -t $NUM_THREADS -m 31 -n $ESTIMATED_GENOME_SIZE -l 31 -n $(($ESTIMATED_GENOME_SIZE*2)) -f `perl -e 'print 1/31/1e5'` mr.fa.in   | grep --text -v '^>' | perl -ane '{$seq=$F[0]; $F[0]=~tr/ACTGactg/TGACtgac/;$revseq=reverse($F[0]); $h{($seq ge $revseq)?$seq:$revseq}=1;}END{$n=0;foreach $k(keys %h){print ">",$n++," length:",length($k),"\n$k\n"}}' > guillaumeKUnitigsAtLeast32bases_all.fasta.tmp && mv guillaumeKUnitigsAtLeast32bases_all.fasta.tmp guillaumeKUnitigsAtLeast32bases_all.31.fasta && \
	    rm -rf work1_mr1 && \
	    createSuperReadsForDirectory.perl --stopAfter joinKUnitigs -minreadsinsuperread 1 -l 31 -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.31.fasta -t $NUM_THREADS -mikedebug work1_mr1 mr.fa.in  1> super1.err 2>&1 && \
            get_super_read_sizes -k work1_mr1/kUnitigLengths.txt -s <(awk '{print $2}' work1_mr1/readPositionsInSuperReads) | sort -S 50% -nrk2 |uniq > work1_mr1/sr_sizes.tmp && \
            reduce_sr `cat work1_mr1/numKUnitigs.txt` work1_mr1/kUnitigLengths.txt 31 work1_mr1/sr_sizes.tmp -o work1_mr1/reduce.tmp 1>reduce2.out 2>&1 && \
            translate_reduced_reads.pl work1_mr1/reduce.tmp < work1_mr1/readPositionsInSuperReads > work1_mr1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt && \
	    find_contained_reads.pl work1_mr1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt $COORDS.1$POSTFIX.fa > containees.txt && \
            trim_by_kunitigs.pl work1_mr1/readPositionsInSuperReads $COORDS.1$POSTFIX.fa work1_mr1/sr_sizes.tmp work1_mr1/kUnitigLengths.txt > $COORDS.1$POSTFIX.trims.txt && \
            trim_mega_reads.pl $COORDS.1$POSTFIX.trims.txt < $COORDS.1$POSTFIX.fa | ufasta extract -v -f containees.txt |make_mr_frg.pl mr 600  > $COORDS.1.frg.tmp && mv  $COORDS.1.frg.tmp  $COORDS.1.frg && \
            trim_mega_reads.pl $COORDS.1$POSTFIX.trims.txt < $COORDS.1$POSTFIX.fa | make_mate_frg.pl > $COORDS.1.mates.frg.tmp && mv $COORDS.1.mates.frg.tmp $COORDS.1.mates.frg && \
            rm -rf $CA work1_mr1 guillaumeKUnitigsAtLeast32bases_all.31.fasta mr.fa.in || error_exit "failed to create mega-reads frg file";
	if  [ ! -s $COORDS.1.frg ];then
	    error_exit "failed to create mega-reads frg file"
	fi
    fi

    TCOVERAGE=20
    if [ $ESTIMATED_GENOME_SIZE -gt 1 ];then
	MR_SIZE=$(stat -c%s -L "$COORDS.1$POSTFIX.fa");
	MCOVERAGE=$(($MR_SIZE/$ESTIMATED_GENOME_SIZE/$PLOIDY+1));
	if [ $MCOVERAGE -le 5 ];then
	    log "Coverage of the mega-reads less than 5 -- using the super reads as well";
	    SR_FRG=$COORDS.sr.frg
	    if [ ! -s $SR_FRG ];then
		awk '{if($0 ~ /^>/) print $0":super-read"; else print $0}' $MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta | fasta2frg.pl sr 200 > $SR_FRG.tmp && mv  $SR_FRG.tmp  $SR_FRG || error_exit "failed to create super-reads frg file";
	    fi
	fi
	COVERAGE=`ls $SR_FRG $COORDS.1.frg $OTHER_FRG 2>/dev/null | xargs stat -c%s | awk '{n+=$1}END{cov=int(n/int('$ESTIMATED_GENOME_SIZE')/int('$PLOIDY')); if(cov<15) cov=15; print cov;}'`;
	TCOVERAGE=$COVERAGE;
    fi

    OVLREFSIZE=`ls $SR_FRG $COORDS.1.frg $OTHER_FRG 2>/dev/null | xargs stat -c%s | perl -ane '$n+=$F[0];END{if(int($n/200)<50000){print "50000";}else{print int($n/200)}}'`

    rm -f .rerun
    rm -f $CA.log

    OVLMIN=`head -n 100000 $SR_FRG $COORDS.1.frg $OTHER_FRG 2>/dev/null | grep -A 1 '^seq:' |grep -v '^seq:' | grep -v '\-\-' | awk 'BEGIN{minlen=100000}{if(length($1)<minlen && length($1)>=64) minlen=length($1);}END{if(minlen>=int("'$OVLMIN_DEFAULT'")) print "'$OVLMIN_DEFAULT'"; else print minlen-1;}'`
    batOptions="-repeatdetect $TCOVERAGE $TCOVERAGE $TCOVERAGE -el $OVLMIN -RS"
    OVL_MER=22

    log "Coverage threshold for splitting unitigs is $TCOVERAGE minimum ovl $OVLMIN"
    let NUM_THREADSd4=$(($NUM_THREADS/4+1))
    if [ $USE_GRID -ge 1 ];then
	OVL_THREADS=4
    else
	OVL_THREADS=2
    fi

    echo "batOptions=$batOptions
useGrid=$USE_GRID
gridEngine=$GRID_ENGINE
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
ovlCorrConcurrency=$NUM_THREADSd4
frgCorrThreads=$NUM_THREADSd4
frgCorrConcurrency=$NUM_THREADSd4
mbtThreads=$NUM_THREADS
ovlThreads=$OVL_THREADS
ovlHashBlockLength=10000000
ovlRefBlockSize=$OVLREFSIZE
ovlConcurrency=$NUM_THREADS
doOverlapBasedTrimming=1
doUnitigSplitting=0
doChimeraDetection=normal
merylThreads=$NUM_THREADS
stoneLevel=0
doExtendClearRanges=0
computeInsertSize=0
maxRepeatLength=12000
ovlErrorRate=0.1
cnsOnGrid=0
cnsConcurrency=$NUM_THREADS
cnsMinFrags=10000
cnsErrorRate=0.1
cnsMaxCoverage=7
cnsReuseUnitigs=1
cgwErrorRate=0.1
cgwMergeMissingThreshold=-1
cgwMergeFilterLevel=1
cgwDemoteRBP=0
cgwPreserveConsensus=1" > runCA.spec


    log "Running assembly"
    if [ ! -e "${CA}/5-consensus/consensus.success" ]; then 
  #need to start from the beginning
  #this is helpful for re-starts
	rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh
	$CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopBefore=scaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1

#sometimes overlap jobs need to be resubmitted -- check of OBT worked
	if [ ! -d "${CA}/1-overlapper" ]; then
	    rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh && \
		$CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopBefore=scaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
	fi

#sometimes overlap jobs need to be resubmitted -- check of OVL worked
	if [ ! -d "${CA}/3-overlapcorrection" ]; then
	    rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh && \
		$CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopBefore=scaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
	fi

#this is a fix for sometimes failing fragment correction
	if [ ! -e "${CA}/4-unitigger/unitigger.success" ]; then
	    rm -f $CA/0-overlaptrim-overlap/overlap.sh $CA/1-overlapper/overlap.sh
	    echo "doFragmentCorrection=0" >> runCA.spec
	    $CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA stopBefore=scaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
	fi
	rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
	$CA_PATH/runCA -s runCA.spec -p genome -d $CA  stopBefore=scaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
    fi

#at athis point we check if the unitig consensus is done
    if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
	error_exit "Assembly stopped or failed, see $CA.log" 
    fi

    if [ ! -e "${CA}/deduplicate.success" ]; then
#here we remove overlaps to the reads in duplicate/redundant unitigs and then re-run the unitigger/consensus
	deduplicate_unitigs.sh $CA_PATH $CA genome $NUM_THREADS $OVL_MER $PLOIDY
    fi

    if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
  #after deduplicate we need to rebuild the unitigs, we rerun CA on deduplicated overlapStore
	$CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopBefore=scaffolder $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
	rm -rf $CA/5-consensus/*.success $CA/5-consensus/consensus.sh
	$CA_PATH/runCA -s runCA.spec -p genome -d $CA  stopBefore=scaffolder cnsConcurrency=$(($NUM_THREADS/2+1)) $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
    fi

    if [ ! -e "${CA}/5-consensus/consensus.success" ]; then
	error_exit "Assembly stopped or failed, see $CA.log"
    fi

#recompute astat if low pacbio coverage
    if [ $MCOVERAGE -le 5 ]; then
	if [ ! -e ${CA}/recompute_astat.success ];then
	    log "Recomputing A-stat"
	    recompute_astat_superreads_CA8.sh genome $CA $PE_AVG_READ_LENGTH $MASURCA_ASSEMBLY_WORK1_PATH/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt  $SR_FRG
	    touch ${CA}/recompute_astat.success
	fi
    fi

#we start from here if the scaffolder has been run or continue here  
    $CA_PATH/runCA -s runCA.spec consensus=pbutgcns -p genome -d $CA  stopBefore=terminator $COORDS.1.frg $COORDS.1.mates.frg $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1
    rm -rf $CA/8-consensus/*.success $CA/8-consensus/consensus.sh
    $CA_PATH/runCA -s runCA.spec -p genome -d $CA  cnsConcurrency=$(($NUM_THREADS/2+1)) $COORDS.1.frg $COORDS.1.mates.frg  $SR_FRG $OTHER_FRG 1>> $CA.log 2>&1 && \
	log "Mega-reads initial assembly complete" || error_exit "Assembly stopped or failed, see $CA.log"

    if [ ! -e "${CA}/10-gapclose/gapclose.success" ] && [ $(stat -c%s ${CA}/9-terminator/genome.scf.fasta) -gt $(stat -c%s ${CA}/9-terminator/genome.ctg.fasta) ] ; then
	log "Closing gaps in scaffolds"
	mkdir -p ${CA}/10-gapclose
	(cd ${CA}/10-gapclose && \
	    $MYPATH/process_scaffold_gaps.pl ../9-terminator/genome.posmap.ctgscf ../9-terminator/genome.posmap.frgctg |sort -k 2 -S 10% > read_scaffold.txt && \
	    $MYPATH/splitScaffoldsAtNs.pl  < ../9-terminator/genome.scf.fasta > genome.scf.split.fa && \
	    grep '^>' --text genome.scf.split.fa | perl -ane '{($rn,$coord)=split(/\./,substr($F[0],1));$h{$rn}.=substr($F[0],1)." ";}END{foreach $r(keys %h){@f=split(/\s+/,$h{$r}); for ($i=0;$i<$#f;$i++){print $f[$i]," ",$f[$i+1],"\n"}}}' > valid_join_pairs.txt && \
	    $MYPATH/ufasta extract -f <(awk '{print $1"\n"$2;}' valid_join_pairs.txt) genome.scf.split.fa > to_join.scf.fa && \
	    $MYPATH/ufasta extract -f <(awk '{print $1;}' read_scaffold.txt) ../../$LONGREADS1 > qrys.all.fa && \
	    $MYPATH/ufasta extract -f <(awk '{if($2 != ps) print $1; ps=$2}' read_scaffold.txt) qrys.all.fa > refs.fa && \
	    perl -ane '{
              $h{$F[0]}=$F[1];
            }END{
              $flag=0;
              open(FILE,"refs.fa");
              while($line=<FILE>){
                if($line=~/^>/){
                  chomp($line);
                  @f=split(/\s+/,$line);
                  if(defined($h{substr($f[0],1)})){
                    $line=<FILE>;
                    chomp($line);
                    $hseq{$h{substr($f[0],1)}}=$line;
                  }
                }
              }
          foreach $name(keys %hseq){
            print ">$name\n$hseq{$name}\n";
          }}' read_scaffold.txt > refs.renamed.fa && \
	  $MYPATH/ufasta extract -v -f <(awk '{if($2 != ps) print $1; ps=$2}' read_scaffold.txt) qrys.all.fa > qrys.fa && \
	  $MYPATH/../CA8/Linux-amd64/bin/blasr qrys.fa  refs.renamed.fa  -nproc $NUM_THREADS -bestn 10 -m 5 2>blasr.err | \
	  sort -k6 -S10% | \
	  $MYPATH/../CA8/Linux-amd64/bin/pbdagcon -j $NUM_THREADS -t 0 -c 1 /dev/stdin  2>pbdagcon.err | \
	  tee join_consensus.fasta | \
	  $MYPATH/nucmer --delta /dev/stdout --batch 10000 -l 17 -c 51 -L 200 -t $NUM_THREADS to_join.scf.fa /dev/stdin | \
	  $MYPATH/show-coords -lcHq /dev/stdin > scf_join.coords && \
	  perl -ane '{($scf1)=split(/\./,$F[-1]);($scf2)=split(/\./,$F[-2]); print if($scf1 eq $scf2);}' scf_join.coords | \
	  $MYPATH/extract_merges_mega-reads.pl join_consensus.fasta  valid_join_pairs.txt > merges.txt && \
	  perl -ane '{if($F[2] eq "F"){$merge="$F[0] $F[3]";}else{$merge="$F[3] $F[0]";} if(not(defined($h{$merge}))|| $h{$merge} > $F[1]+$F[4]){$hl{$merge}=join(" ",@F);$h{$merge}=$F[1]+$F[4];}}END{foreach $k(keys %hl){print $hl{$k},"\n"}}' merges.txt > merges.best.txt && \
	  cat <($MYPATH/ufasta extract -v -f <(awk '{print $1"\n"$2;}' valid_join_pairs.txt) genome.scf.split.fa) <($MYPATH/merge_mega-reads.pl < merges.best.txt | $MYPATH/create_merged_mega-reads.pl to_join.scf.fa  merges.best.txt) | recover_scaffolds.pl > genome.scf.joined.fa.tmp && mv genome.scf.joined.fa.tmp genome.scf.joined.fa && \
          $MYPATH/samba.sh -r genome.scf.joined.fa -q $LONGREADS -t $NUM_THREADS && \
          mv genome.scf.joined.fa.scaffolds.fa genome.scf.fasta && touch gapclose.success
	)
    fi
fi #if flye
#scaffolding with the long reads

