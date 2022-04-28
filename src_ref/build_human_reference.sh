#!/bin/bash
#chm13.draft_v1.1.mat.fasta.w_noise  chm13.draft_v1.1.pat.fasta.w_noise  chm13.draft_v1.1.withY.fasta  chm13.draft_v1.1.withY.fasta.w_noise
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
PAT_ASM="paternal.fa"
MAT_ASM="maternal.fa"
HAP_ASM="haploid.fa"
ALT_ASM="alternative.fa"
READS=""
NUM_THREADS=2
CHM13PATH="/"
CHM13MALE=chm13.draft_v1.1.withY.fasta
CHM13=chm13.draft_v1.1.fasta
CHM13PAT=chm13.draft_v1.1.pat.fasta
CONTAMINANTS=contaminants.fa
GENDER="none"

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


function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage () {
echo "Usage:"
echo "-t <number of threads, default:1>"
echo "-p <paternal genome fasta>"
echo "-m <maternal genome fasta>"
echo "-b <haplotype merged assembly fasta>"
echo "MANDATORY: either -b hap_merged.fasta or (-m maternal.fasta and -p paternal.fasta) must be supplied"
echo "-g <MANDATORY gender: male or female>"
echo "-c <MANDATORY: path to chm13 reference files downloaded from ftp://ftp.ccb.jhu.edu/pub/alekseyz/chm13.tgz"
echo "-r <MANDATORY: reads fastq or fasta>"
echo "-v verbose mode"
echo "-h this usage message"
}

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            NUM_THREADS="$2"
            shift
            ;;
        -c|--chm13)
            CHM13PATH="$2"
            shift
            ;;
        -p|--paternal)
            PAT_ASM="$2"
            shift
            ;;
        -m|--maternal)
            MAT_ASM="$2"
            shift
            ;;
        -g|--gender)
            GENDER="$2"
            shift
            ;;
        -b|--both-haplotypes)
            HAP_ASM="$2"
            shift
            ;;
        -r|--reads)
            READS="$2"
            shift
            ;;
        -a|--alt-asm)
            ALT_ASM="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [ ! -s $CHM13PATH/$CHM13 ] || [ ! -s $CHM13PATH/$CHM13MALE ] || [ ! -s $CHM13PATH/$CHM13PAT.w_noise ] || [ ! -s $CHM13PATH/$CONTAMINANTS ] || [ ! -s $CHM13PATH/chrM.fa ];then
  error_exit "reference genome files not found, please download them from ftp://ftp.ccb.jhu.edu/pub/alekseyz/chm13.tgz, unpack (tar -xvzf chm13.tgz), and specify the PATH to them with -c /path_to/"
fi

if [ ! $GENDER = "male" ] && [ ! $GENDER = "male" ];then
  error_exit "please specify gender with -g male|female"
fi

if [ -s $READS ];then
  #make all paths absolute
  READS_FN=`basename $READS`
  DPATH="`dirname \"$READS\"`"
  DPATH="`( cd \"$DPATH\" && pwd )`"
  READS=$DPATH/$READS_FN
else
  error_exit "Reads file $READS is not provided or does not exist.  Please provide it with -r <reads.fastq>.  this argumant is mandatory"
fi

if [ ! -e build.chromosomes.success ];then
  if [ -s $HAP_ASM ];then
      log "building chromosomes from merged haplotype assembly" && \
      mkdir -p haploid && \
      HAP_ASM_FN=`basename $HAP_ASM` && \
      DPATH="`dirname \"$HAP_ASM\"`" && \
      DPATH="`( cd \"$DPATH\" && pwd )`" && \
      HAP_ASM=$DPATH/$HAP_ASM_FN && \
      (cd haploid && \
      if [ $GENDER = "male" ];then
        REF=$CHM13MALE
      else
        REF=$CHM13
      fi && \
      chromosome_scaffolder.sh -c 50000 -i 99 -m 250000 -r $CHM13PATH/$REF.w_noise -q $HAP_ASM -s $READS -t $NUM_THREADS -n 0 && \
      close_scaffold_gaps.sh -r $REF.w_noise.$HAP_ASM_FN.split.reconciled.fa -q $READS -t $NUM_THREADS) && \
      if [ -s haploid/$REF.w_noise.$HAP_ASM_FN.split.reconciled.fa.split.joined.fa ];then
        cat haploid/$REF.w_noise.$HAP_ASM_FN.split.reconciled.fa.split.joined.fa  > chromosomes.fa.tmp && \
        mv chromosomes.fa.tmp chromosomes.fa
      else
        error_exit "error building chromosomes: check haploid/$REF.w_noise.$HAP_ASM_FN.split.reconciled.fa.split.joined.fa"
      fi      
  elif [ -s $MAT_ASM ] && [ -s $PAT_ASM ];then
    log "building chromosomes from haplotype resolved assembly" && \
    mkdir -p maternal && \
    MAT_ASM_FN=`basename $MAT_ASM` && \
    DPATH="`dirname \"$MAT_ASM\"`" && \
    DPATH="`( cd \"$DPATH\" && pwd )`" && \
    MAT_ASM=$DPATH/$MAT_ASM_FN && \
    (cd maternal && \
    chromosome_scaffolder.sh -c 50000 -i 99 -m 250000 -r $CHM13PATH/$CHM13.w_noise -q $MAT_ASM -s $READS -t $NUM_THREADS -n 0 && \
    close_scaffold_gaps.sh -r $CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa -q $READS -t $NUM_THREADS) && \
    if [ $GENDER = "male" ];then

      mkdir -p paternal
      PAT_ASM_FN=`basename $PAT_ASM` && \
      DPATH="`dirname \"$PAT_ASM\"`" && \
      DPATH="`( cd \"$DPATH\" && pwd )`" && \
      PAT_ASM=$DPATH/$PAT_ASM_FN && \
      (cd paternal && \
      chromosome_scaffolder.sh -c 50000 -i 99 -m 250000 -r $CHM13PATH/$CHM13PAT.w_noise -q $PAT_ASM -s $READS -t $NUM_THREADS -n 0 && \
      close_scaffold_gaps.sh -r $CHM13PAT.w_noise.$PAT_ASM_FN.split.reconciled.fa -q $READS -t $NUM_THREADS) && \
      if [ -s maternal/$CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa.split.joined.fa ]  && [  -s paternal/$CHM13PAT.w_noise.$PAT_ASM_FN.split.reconciled.fa.split.joined.fa ];then
        cat maternal/$CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa.split.joined.fa  <(ufasta extract -n chrY paternal/$CHM13PAT.w_noise.$PAT_ASM_FN.split.reconciled.fa.split.joined.fa ) > chromosomes.fa.tmp && \
        mv chromosomes.fa.tmp chromosomes.fa
      else
        error_exit "error building chromosomes: check maternal/$CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa.split.joined.fa and paternal/$CHM13PAT.w_noise.$PAT_ASM_FN.split.reconciled.fa.split.joined.fa"
      fi
    else #if female
      if [ -s maternal/$CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa.split.joined.fa ];then
        cat maternal/$CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa.split.joined.fa  > chromosomes.fa.tmp && \
        mv chromosomes.fa.tmp chromosomes.fa 
      else
        error_exit "error building chromosomes: check maternal/$CHM13.w_noise.$MAT_ASM_FN.split.reconciled.fa.split.joined.fa"
      fi
    fi #if GENDER
  else
    usage
    error_exit "Either maternal and paternal haplotype or haplotype-merged scaffolds must be supplied)"
  fi

  if [ -s $ALT_ASM ];then
    log "Closing gaps with an alternative assembly"
    ALT_ASM_FN=`basename $ALT_ASM` && \
    DPATH="`dirname \"$ALT_ASM\"`" && \
    DPATH="`( cd \"$DPATH\" && pwd )`" && \
    ALT_ASM=$DPATH/$ALT_ASM_FN && \
    mkdir -p alternative && \
    (cd alternative && \
    close_scaffold_gaps.sh -d asm -r ../chromosomes.fa -q $ALT_ASM -t $NUM_THREADS -m 2500 -o 10000 1>out.err 2>&1) && \
    mv chromosomes.fa chromosomes.fa.bak && cp alternative/chromosomes.fa.split.joined.fa chromosomes.fa.tmp && \
    mv chromosomes.fa.tmp chromosomes.fa
  fi
  touch build.chromosomes.success && rm -f build.gapless.success
fi

if [ ! -e  build.gapless.success ];then
  log "building gapless chromosomes"
  echo "#!/bin/bash" > run.sh && \
  for c in $(seq 1 22;echo " X");do
    echo "bash $MYPATH/final_polish.sh $c $CHM13PATH/$CHM13 $CHM13PATH/$CHM13.w_noise chromosomes.fa &" >>run.sh
  done
  echo "wait" >> run.sh
  bash ./run.sh && \
  log "Building chrM" && \
  ufasta extract -f <(ufasta sizes -H chromosomes.fa | grep chrM |awk '{print $1}') chromosomes.fa > chrM.draft.fa && \
  $MYPATH/extract_chrM.sh $CHM13PATH/chrM.fa chrM.draft.fa && \
  if [ $GENDER = "male" ];then
    log "Building Y" && \
    mkdir -p Y.dir && \
    if [ ! -s Y.dir/Y.all.polished.fa ];then
      (cd Y.dir && \
      ufasta extract -n chrY $CHM13PATH/$CHM13MALE > refY.fa.tmp && mv refY.fa.tmp refY.fa && \
      ufasta extract -n chrY ../chromosomes.fa > asmY.fa.tmp && mv asmY.fa.tmp asmY.fa && \
      ufasta one refY.fa |perl -ane 'BEGIN{$keep_on_ends=100000;}{if($F[0] =~ /^>/){print}else{$sub=length($F[0])-2*$keep_on_ends;if($sub>10000){substr($F[0],$keep_on_ends,$sub,"N"x$sub); print $F[0],"\n";}else{print}}}' | \
      splitScaffoldsAtNs.pl > refY.BEG_END.fa && \
      cat <(perl -e '{print ">chrY\n"}') \
      <(grep -v '^>' refY.BEG_END.fa |head -n 1) \
      <(perl -e '{print "N"x1000}') \
      <(grep -v '^>' asmY.fa ) \
      <(perl -e '{print "N"x1000}') \
      <(grep -v '^>' refY.BEG_END.fa |tail -n 1) > asmY.padded.fa.tmp && mv asmY.padded.fa.tmp asmY.padded.fa && \
      splitScaffoldsAtNs.pl  < asmY.padded.fa | ufasta one > asmY.split.fa.tmp && mv asmY.split.fa.tmp asmY.split.fa && \
      grep '^>' --text asmY.split.fa | perl -ane '{($rn,$coord)=split(/\./,substr($F[0],1));$h{$rn}.=substr($F[0],1)." ";}END{foreach $r(keys %h){@f=split(/\s+/,$h{$r}); for ($i=0;$i<$#f;$i++){print $f[$i]," ",$f[$i+1],"\n"}}}' > asmY.split.fa.valid_join_pairs.tmp && mv asmY.split.fa.valid_join_pairs.tmp asmY.split.fa.valid_join_pairs && \
      nucmer -p aln -l 21 -c 200 -t $NUM_THREADS <(ufasta one asmY.split.fa |perl -ane 'BEGIN{$keep_on_ends=1000000;}{if($F[0] =~ /^>/){print}else{$sub=length($F[0])-2*$keep_on_ends;if($sub>10000){substr($F[0],$keep_on_ends,$sub,"N"x$sub); print $F[0],"\n";}else{print}}}') refY.fa && \
      delta-filter -1 aln.delta > aln.1delta.tmp && mv aln.1delta.tmp aln.1delta && \
      show-coords -I 99.5 -lcHq aln.1delta | extract_merges.pl refY.fa 2500 1000000 asm asmY.split.fa.valid_join_pairs > Y.links.tmp && mv Y.links.tmp Y.links && \
      merge_contigs.pl asmY.split.fa < Y.links | create_merged_sequences.pl asmY.split.fa  Y.links > Y.scaffolds.fa.tmp && mv Y.scaffolds.fa.tmp Y.scaffolds.fa && \
      recover_scaffolds.pl < Y.scaffolds.fa |ufasta format > Y.joined.fa.tmp && mv Y.joined.fa.tmp Y.all.polished.fa )
    fi
  fi
  cat *.dir/*.all.polished.fa chrM.new.fa > chromosomes.gapless.fa.tmp && \
  rm chrM.draft.fa chrM.new.fa && \
  mv chromosomes.gapless.fa.tmp chromosomes.gapless.fa &&  \
  touch build.gapless.success && rm -f build.add_singletons.success
fi #if success

if [ ! -e build.add_singletons.success ];then
  log "adding unplaced and screening for contaminants"
  ufasta extract -f <(ufasta sizes -H chromosomes.fa | grep -v chr |awk '{print $1}') chromosomes.fa > chromosomes.unplaced.fa.tmp && \
  mv chromosomes.unplaced.fa.tmp chromosomes.unplaced.fa && \
  nucmer -t $NUM_THREADS -p contaminant chromosomes.unplaced.fa $CHM13PATH/$CONTAMINANTS && \
  nucmer -t $NUM_THREADS -p unplaced chromosomes.unplaced.fa chromosomes.gapless.fa && \
  ufasta extract -v -f <(cat <(delta-filter -r -i 95 unplaced.delta |show-coords -lcHr /dev/stdin | perl -ane '$pct{$F[-2]}+=$F[14];END{foreach $k (keys %pct){print "$k $pct{$k}\n";}}' | awk '{if($2>=75) print $1}') <(delta-filter -r -i 95 contaminant.delta |show-coords -lcHr /dev/stdin | perl -ane '$pct{$F[-2]}+=$F[14];END{foreach $k (keys %pct){print "$k $pct{$k}\n";}}' | awk '{if($2>=50) print $1}')) chromosomes.unplaced.fa > chromosomes.to_add.fa.tmp && \
  cat chromosomes.gapless.fa chromosomes.to_add.fa.tmp > chromosomes.gapless.w_unplaced.fa.tmp && \
  mv chromosomes.gapless.w_unplaced.fa.tmp chromosomes.gapless.w_unplaced.fa && rm chromosomes.to_add.fa.tmp && \
  touch build.add_singletons.success
fi

if [ ! -e build.polish.success ];then
  log "polishing sequence"
  /ccb/salz4-2/alekseyz/mer_polish/build/bin/jasper.sh -t 24 -b 800000000 -a chromosomes.gapless.w_unplaced.fa -r $READS -k 25 && \
  touch build.polish.success
fi

if [ -e build.polish.success ];then
  log "output assembly ready for annotation is in chromosomes.gapless.w_unplaced.fa.fixed.fa"
fi


