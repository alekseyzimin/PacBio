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
USE_GRID=0
PACBIO=""
NANOPORE=""
NANOPORE_RNA=0
ONEPASS=0
OVLMIN_DEFAULT=250
FLYE=0
POSTFIX=""
#MIN_PROPORTION="0.15"
#MIN_RADIUS="15"
MIN_PROPORTION="0.25"
MIN_RADIUS="400"
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
        -Pb|--pbatch-size)
            PBATCH_SIZE="$2"
            shift
            ;;
	-q|--queue)
	    QUEUE="$2"
	    shift
	    ;;
        -1|--onepass)
            ONEPASS=1;
            ;;
	-o|--other_frg)
	    OTHER_FRG="$2"
	    shift
	    ;;
        -F|--Flye)
            FLYE=1;
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
if [ $FLYE -gt 0 ];then
    log "Using Flye from $CA_PATH"
    if [ ! -e $CA_PATH/flye ];then
	error_exit "flye not found at $CA_PATH!";
    fi
#    POSTFIX=".flye"
#    MIN_PROPORTION="0.25"
#    MIN_RADIUS="400"
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
COORDS=mr.$KMER.$MER.$B.$d
CA=CA.${COORDS}
echo $CA > CA_dir.txt
echo $PLOIDY > PLOIDY.txt

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


zcat -f $LONGREADS | $MYPATH/fastqTofasta.pl > $LONGREADS1.tmp && mv $LONGREADS1.tmp $LONGREADS1
LONGREADS1="long_reads.fa";
SUPERREADS=superReadSequences.named.fasta
KUNITIGS=$MASURCA_ASSEMBLY_WORK1_PATH/../guillaumeKUnitigsAtLeast32bases_all.fasta
perl -ane 'push(@names,$F[0]);
END{
  open(FILE,"'$MASURCA_ASSEMBLY_WORK1_PATH/superReadSequences.fasta.all'");
  while($line=<FILE>){
    if($line=~/^>/){
      chomp($line);
      print ">",$names[substr($line,1)],"\n";
    }else{
      print $line;
    }
  }
}' < $MASURCA_ASSEMBLY_WORK1_PATH/superReadNames.txt > $SUPERREADS.tmp && mv $SUPERREADS.tmp $SUPERREADS || error_exit "failed to create named super-reads file";

if [ ! -s $SUPERREADS ];then
    error_exit "Error creating named super-reads file ";
fi

if [ ! -s $KUNITIGS ];then
    error_exit "K-unitigs file $KUNITIGS not found; failed to create named super-reads file";
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
    log "Mega-reads pass 1"

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
		    echo "$MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SGE_TASK_ID -o mr.batch\$SGE_TASK_ID.tmp && mv mr.batch\$SGE_TASK_ID.tmp mr.batch\$SGE_TASK_ID.txt && touch mr.batch\$SGE_TASK_ID.success" >> create_mega_reads.sh && \
		    echo "else" >> create_mega_reads.sh && \
		    echo "echo \"job \$SGE_TASK_ID previously completed successfully\"" >> create_mega_reads.sh && \
		    echo "fi"  >> create_mega_reads.sh && chmod 0755 create_mega_reads.sh
            else
		echo "#!/bin/sh" > create_mega_reads.sh && \
		    echo "if [ ! -e mr.batch\$SLURM_ARRAY_TASK_ID.success ];then" >> create_mega_reads.sh && \
		    echo "$MYPATH/create_mega_reads -s $JF_SIZE -m $MER --psa-min 12  --stretch-cap 10000 -k $KMER -u ../$KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r ../$SUPERREADS  -p lr.batch\$SLURM_ARRAY_TASK_ID -o mr.batch\$SLURM_ARRAY_TASK_ID.tmp && mv mr.batch\$SLURM_ARRAY_TASK_ID.tmp mr.batch\$SLURM_ARRAY_TASK_ID.txt && touch mr.batch\$SLURM_ARRAY_TASK_ID.success" >> create_mega_reads.sh && \
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
            numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p <(ufasta extract -v -f <(head -n -1 $COORDS.txt.headers | awk '{print substr($2,2)}') $LONGREADS1)  -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt 
          else
            create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p <(ufasta extract -v -f <(head -n -1 $COORDS.txt.headers | awk '{print substr($2,2)}') $LONGREADS1)  -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt
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
            numactl --interleave=all create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed"
          else
            create_mega_reads -s $JF_SIZE -m $MER --psa-min 13  --stretch-cap 10000 -k $KMER -u $KUNITIGS -t $NUM_THREADS -B $B --max-count 5000 -d $d  -r $SUPERREADS  -p $LONGREADS1 -o $COORDS.txt.tmp 1>create_mega-reads.err 2>&1 && mv $COORDS.txt.tmp $COORDS.txt || error_exit "mega-reads pass 1 failed"
          fi
        fi #failed run
    fi #single computer
    touch .rerun
    if  [ ! -s $COORDS.txt ];then
	error_exit "mega-reads pass 1 failed"
    fi
fi

touch $COORDS.mr.txt
if [ ! -s $COORDS.single.txt ] || [ -e .rerun ];then
  grep --text '^>' $COORDS.txt | awk  '{print substr($1,2)}' > $COORDS.single.txt
fi

if [ ! -s ${COORDS}.all.txt ] || [ -e .rerun ];then
    log "Refining alignments"
    NUM_LONGREADS_READS_PER_BATCH=`grep --text '^>'  $LONGREADS1 | wc -l | awk '{bs=int($1/1024);if(bs<1000){bs=1000};if(bs>100000){bs=100000};}END{print bs}'` 
    cat <(ufasta extract -f $COORDS.single.txt $COORDS.txt) <(ufasta extract -v -f $COORDS.single.txt $COORDS.mr.txt)| awk '{if($0~/^>/){pb=substr($1,2);print $0} else { print $3" "$4" "$5" "$6" "$10" "pb" "$11" "$9}}' | add_pb_seq.pl $LONGREADS1 | split_matches_file.pl $NUM_LONGREADS_READS_PER_BATCH .matches && ls .matches.* | xargs -P $NUM_THREADS -I % refine.sh $COORDS % $KMER && cat $COORDS.matches*.all.txt.tmp > $COORDS.all.txt && rm .matches.* && rm $COORDS.matches*.all.txt.tmp 
    touch .rerun
fi

#for now we just exit here
join_mega_reads_trim.onepass.ref.pl <$COORDS.all.txt 1>$COORDS.transcripts.fa 2>$COORDS.transcripts.err




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
	ufasta split -i <(add_pb_seq.pl $LONGREADS1 < ${COORDS}.all.txt) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.1 2>$COORDS.1$POSTFIX.to_join.fa.tmp.1) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.2 2>$COORDS.1$POSTFIX.to_join.fa.tmp.2) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.3 2>$COORDS.1$POSTFIX.to_join.fa.tmp.3) \
            >($MYPATH/join_mega_reads_trim.onepass.nomatch.pl ${COORDS}.1$POSTFIX.allowed  $MAX_GAP 1>$COORDS.1$POSTFIX.fa.tmp.4 2>$COORDS.1$POSTFIX.to_join.fa.tmp.4) && \
            cat $COORDS.1$POSTFIX.fa.tmp.{1,2,3,4} | awk '{if($1!~/^>/ && $1~/>/){split($1,a,">");print a[1];if(a[2]!="") print ">"a[2];}else{print $0}}' > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.unjoined.fa && \
            cat $COORDS.1$POSTFIX.to_join.fa.tmp.{1,2,3,4} > $COORDS.1$POSTFIX.to_join.fa.tmp && mv $COORDS.1$POSTFIX.to_join.fa.tmp $COORDS.1$POSTFIX.to_join.fa 
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
                    log "Using SGE grid queue $QUEUE"
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
                  echo "sort -k6 -S2% | $MYPATH/../CA8/Linux-amd64/bin/pbdagcon -j $NUM_THREADS -t 0 -c 1 /dev/stdin  2>pbdagcon.err | awk -F 'N' '{if(\$1 == \"\") print \"ACGT\"; else print \$1}' > join_consensus.\$TASK_ID.fasta && \\" >> ./do_consensus.sh
                  echo "$MYPATH/nucmer --delta /dev/stdout --batch 10000 -l 17 -c 51 -L 200 -t $NUM_THREADS to_join.\$TASK_ID.fa join_consensus.\$TASK_ID.fasta 2>/dev/null | \\" >> ./do_consensus.sh
                else
                  echo "$MYPATH/../CA8/Linux-amd64/bin/blasr to_blasr.\$TASK_ID.fa   ref.\$TASK_ID.fa  -minMatch 15 -nproc 16 -bestn 10 -m 5 2>blasr.err | \\" >> ./do_consensus.sh 
                  echo "sort -k6 -S2% | $MYPATH/../CA8/Linux-amd64/bin/pbdagcon -j 8 -t 0 -c 1 /dev/stdin  2>pbdagcon.err | awk -F 'N' '{if(\$1 == \"\") print \"ACGT\"; else print \$1}' > join_consensus.\$TASK_ID.fasta && \\" >> ./do_consensus.sh
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
                    qsub -q $QUEUE -cwd -j y -sync y -N "join_mega_reads"  -t 1-$TOJOIN_BATCHES do_consensus.sh 1> cqsub2.out 2>&1 || error_exit "join consensus failed on the grid"  
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
            cat $COORDS.1$POSTFIX.joined.fa $COORDS.1$POSTFIX.unjoined.fa  > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.fa
        else
            cat $COORDS.1$POSTFIX.unjoined.fa $COORDS.1$POSTFIX.to_join.fa  > $COORDS.1$POSTFIX.fa.tmp && mv $COORDS.1$POSTFIX.fa.tmp $COORDS.1$POSTFIX.fa
        fi
    fi
    touch .rerun
    if  [ ! -s $COORDS.1$POSTFIX.fa ];then
        error_exit "Gap consensus failed"
    fi
fi

