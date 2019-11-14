#!/bin/bash
#this code aims at reconciling the hybrid contigs and the chromosomes of the previously produces assembly
#arguments are: reference chromosomes, hybrid contigs, hybrid posmap (frgctg), filtered delta-file of alignments of ref to hyb
#MUST HAVE blasr on the PATH
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"

export PATH=$MYPATH:$PATH;
set -o pipefail
set -e
NUM_THREADS=1
#default minimum alignment identity
IDENTITY=97
#parameter for merging alignments
MERGE=100000
MERGE_SEQ=0

#low coverage threshold for breaking
COV_THRESH=3
REP_COV_THRESH=30
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


#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            NUM_THREADS="$2"
            shift
            ;;
        -s|--sequenced_reads)
            READS="$2"
            shift
            ;;
        -q|--query)
            QRY="$2"
            shift
            ;;
        -i|--identity)
            IDENTITY="$2"
            shift
            ;;
        -cl|--low_coverage_threshold)
            COV_THRESH="$2"
            shift
            ;;
        -ch|--repeat_threshold)
            REP_COV_THRESH="$2"
            shift
            ;;
        -m|--merge-slack)
            MERGE="$2"
            shift
            ;;
        -M|--merge-sequences)
            MERGE_SEQ=1
            ;;
        -r|--reference)
            REF="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: chromosome_scaffolder.sh -r <reference genome> -q <assembly to be scaffolded with the reference> -t <number of threads> -i <minimum sequence similarity percentage> -m <merge polishing sequence alignments slack (advanced)> -v <verbose> -s <reads to align to the assembly to check for misassemblies> -cl <coverage threshold for splitting at misassemblies, default 3> -ch <repeat coverage threshold for splitting at misassemblies, default 30>" 
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

REF_CHR=`basename $REF`
HYB_CTG=`basename $QRY`.split
HYB_POS=$HYB_CTG.posmap
PREFIX=$REF_CHR.$HYB_CTG
let IDENTITY=$IDENTITY-1

if [ ! -e $PREFIX.split.success ];then
  log "Splitting query scaffolds into contigs"
  rm -f $PREFIX.blasr.success
  $MYPATH/splitFileAtNs $QRY 1 > $HYB_CTG   && touch $PREFIX.split.success
fi

if [ ! -e $PREFIX.blasr.success ];then
  log "Mapping reads to query contigs"
  rm -f $PREFIX.coverage.success
  if [[ $READS = *.fa ]] || [[ $READS = *.fasta ]] || [[ $READS = *.fastq ]];then
  $MYPATH/../CA8/Linux-amd64/bin/blasr -nproc $NUM_THREADS -bestn 1 $READS $HYB_CTG | awk '{if(($11-$10)/$12>0.75){if($4==0) print $1" "substr($2,4)" "$7" "$8" f"; else print  $1" "substr($2,4)" "$9-$8" "$9-$7" r"}}' $1 | sort -nk2 -k3n -S 10% > $HYB_POS && touch $PREFIX.blasr.success
  else
  error_exit "Wrong type/extension for the $READS file, must be .fa, .fasta or .fastq"
  fi
fi

if [ ! -e $PREFIX.noise.success ];then
  log "adding noise to reference to align to duplicated regions"
  rm -f $PREFIX.align1.success
  rm -f $PREFIX.align2.success
  $MYPATH/introduce_errors_fasta_file.pl $REF 0.01 1 $REF | $MYPATH/fix_consensus_from_vcf.pl $REF > $REF.w_noise && touch $PREFIX.noise.success
fi 

if [ ! -e $PREFIX.align1.success ];then
  log "Aligning query contigs to reference scaffolds"
  rm -f $PREFIX.filter1.success
  #$MYPATH/nucmer -t $NUM_THREADS -p $REF_CHR.$HYB_CTG -c 200 $REF.w_noise $HYB_CTG && touch $PREFIX.align1.success
  $MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -k 21 -a -Q $REF $HYB_CTG 2>minimap2.err | $MYPATH/samToDelta > $REF_CHR.$HYB_CTG.delta && touch $PREFIX.align1.success
fi

if [ ! -e $PREFIX.filter1.success ];then
  log "Filtering the alignments" 
  rm -f $PREFIX.merge1.success
  $MYPATH/delta-filter -1 -i $IDENTITY -o 20 $REF_CHR.$HYB_CTG.delta >$REF_CHR.$HYB_CTG.1.delta && touch $PREFIX.filter1.success
fi

#compute coverage from the posmap file
if [ ! -e $PREFIX.coverage.success ];then
  log "Computing read coverage for query contigs" 
  rm -f $PREFIX.break.success
  awk '{print $1" "$2" "$3"\n"$1" "$2" "$4}' $HYB_POS | grep -v F |grep -v R | sort -nk2 -k3n -S 10% | $MYPATH/compute_coverage.pl > $HYB_POS.coverage && touch $PREFIX.coverage.success
fi

if [ ! -e $PREFIX.gaps.success ];then
  log "Computing gap coordinates in the reference"
  rm -f $PREFIX.scaffold.success
  $MYPATH/splitFileAtNs $REF 1 > $REF_CHR.split.fasta
  perl -ane '{$h{substr($F[1],3)}=$F[0]}END{while($line=<STDIN>){chomp($line);@f=split(/\s+/,$line);print "$f[0] $h{$f[1]} ",$f[2]+1," $f[3] $f[4]\n";}}' scaffNameTranslations.txt < genome.posmap.ctgscf | awk 'BEGIN{pg=0}{print $2" "pg" "$3;pg=$4}' > $PREFIX.gap_coordinates.txt && touch $PREFIX.gaps.success
fi

if [ ! -e $PREFIX.merge1.success ];then
  log "Merging alignments"
  rm -f $PREFIX.break.success
  show-coords -lcHr $REF_CHR.$HYB_CTG.1.delta | \
  $MYPATH/merge_matches_and_tile_coords_file.pl $MERGE | \
  awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
  awk '{if($16>5 || $7>5000 ) print $0}' > $REF_CHR.$HYB_CTG.1.coords  && touch $PREFIX.merge1.success
fi

if [ ! -e $PREFIX.break.success ];then
  log "Splitting query contigs at suspect locations"
  rm -f $PREFIX.align2.success
  awk '{if($4<$5) print $4" "$5" "($4+$5)/2" "$NF" "$13; else print $5" "$4" "($4+$5)/2" "$NF" "$13;}' $REF_CHR.$HYB_CTG.1.coords| \
  sort -k4 -k1n | \
  uniq -D -f 3 | \
  awk '{if($NF != prev){offset=$2;prev=$NF;print $0}else if($2>offset){print $0;offset=$2;}else{print "contained "$0}}' > $REF_CHR.$HYB_CTG.1.coords.split_contain && \
  grep -v ^contained $REF_CHR.$HYB_CTG.1.coords.split_contain | uniq -D -f 3 | awk '{if($4==ctg){if($1>5000 && $1<$NF-5000) print "alnbreak "substr($4,4)" "$1" 0"}else{ctg=$4;if($2>5000 && $2<$NF-5000) print "alnbreak "substr($4,4)" "$2" 0"}}' > $REF_CHR.$HYB_CTG.1.coords.breaks && \
  cat $REF_CHR.$HYB_CTG.1.coords.breaks $HYB_POS.coverage  | sort -nk2 -k3n -S 10% > $HYB_POS.coverage.w_breaks && \
  grep -C 50 break $HYB_POS.coverage.w_breaks  | $MYPATH/evaluate_splits.pl <(ufasta sizes -H $HYB_CTG | awk '{print substr($1,4)" "$2}') | sort -nk3 -S 10% >  $HYB_POS.coverage.w_breaks.validated && \
  $MYPATH/break_contigs.pl <(grep -v "end" $HYB_POS.coverage.w_breaks.validated |awk '{if($4<=int("'$COV_THRESH'") || ($1="repeat" && $4>=int("'$REP_COV_THRESH'"))) print $0}') < $HYB_CTG > $HYB_CTG.broken && touch $PREFIX.break.success
fi

#now we re-align the broken contigs to the reference
if [ ! -e $PREFIX.align2.success ];then
  log "Re-aligning contigs after splitting"
  rm -f $PREFIX.filter2.success
  #$MYPATH/nucmer --batch 100000000 -t $NUM_THREADS -p $REF_CHR.$HYB_CTG.broken -c 200  $REF.w_noise $HYB_CTG.broken && touch $PREFIX.align2.success
  $MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -k 21 -a -Q $REF $HYB_CTG.broken 2>minimap2.err | $MYPATH/samToDelta > $REF_CHR.$HYB_CTG.broken.delta && touch $PREFIX.align2.success
fi

if [ ! -e $PREFIX.filter2.success ];then
  log "Filtering the alignments"
  rm -f $PREFIX.scaffold.success
  $MYPATH/delta-filter -r -o 99 -i $IDENTITY $REF_CHR.$HYB_CTG.broken.delta | $MYPATH/delta-filter -q /dev/stdin > $REF_CHR.$HYB_CTG.broken.1.delta && touch $PREFIX.filter2.success
fi

#now we merge/rebuild chromosomes
if [ ! -e $PREFIX.scaffold.success ];then
  rm -f $PREFIX.place_extra.success
  log "Final scaffolding"
  if [ $MERGE_SEQ -gt 0 ];then
    $MYPATH/show-coords -lcHr $REF_CHR.$HYB_CTG.broken.1.delta | \
    $MYPATH/merge_matches_and_tile_coords_file.pl $MERGE | \
    awk '{if($16>25 || $7>2000 ) print $0}' |\
    awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
    $MYPATH/extract_single_best_match_coords_file.pl  |\
    $MYPATH/fill_unaligned_gaps.pl $REF 2>$PREFIX.fillseq.fa | perl -ane '{$dir="f"; $dir="r" if($F[3]>$F[4]);print "$F[-2] $F[-1] 1 $F[12] $dir 100 100 $F[12]\n"}' |\
    $MYPATH/output_reconciled_scaffolds.pl <(cat $PREFIX.fillseq.fa $HYB_CTG.broken) |\
    ufasta format > $REF_CHR.$HYB_CTG.reconciled.fa && touch $PREFIX.scaffold.success && touch $PREFIX.place_extra.success
  else
    $MYPATH/show-coords -lcHr $REF_CHR.$HYB_CTG.broken.1.delta | \
    $MYPATH/merge_matches_and_tile_coords_file.pl $MERGE | \
    awk '{if($16>25 || $7>2000 ) print $0}' |\
    awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
    $MYPATH/extract_single_best_match_coords_file.pl  |\
    $MYPATH/reconcile_matches.pl $PREFIX.gap_coordinates.txt  > $PREFIX.reconciled.txt && touch $PREFIX.scaffold.success 
  fi
fi

#we attempt to place extra contigs that were left unmapped due to repeats
if [ ! -e $PREFIX.place_extra.success ];then
  log "Placing extra repeats"
  $MYPATH/ufasta extract -v -f <(awk '{print $2}' $PREFIX.reconciled.txt) $HYB_CTG.broken > $PREFIX.unplaced.fa && \
  $MYPATH/../Flye/bin/flye-minimap2 -t $NUM_THREADS -k 21 -a -Q $PREFIX.unplaced.fa $REF 2>minimap2.err | $MYPATH/samToDelta > $PREFIX.map_unplaced.delta 
  #$MYPATH/nucmer -t $NUM_THREADS  --maxmatch --batch 10000000 -b 100 -l 100 -c 100 -p $PREFIX.map_unplaced $PREFIX.unplaced.fa $REF 
  if [ -s $PREFIX.map_unplaced.delta ];then
    cat $REF_CHR.$HYB_CTG.broken.1.delta <(awk '{if($0 ~ /^>/){print ">"$2" "substr($1,2)" "$4" "$3}else if(NF==7){if($3<$4){print $3" "$4" "$1" "$2" "$6" "$5" "$7}else{print $4" "$3" "$2" "$1" "$6" "$5" "$7}}else{print $0}}' $PREFIX.map_unplaced.delta | $MYPATH/delta-filter -r -o 99 -i $IDENTITY /dev/stdin | $MYPATH/delta-filter -q /dev/stdin | tail -n +3 ) |\
    $MYPATH/show-coords -lcHr /dev/stdin | \
    $MYPATH/merge_matches_and_tile_coords_file.pl $MERGE | \
    awk '{if($16>25 || $7>2000 ) print $0}' |\
    awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
    $MYPATH/extract_single_best_match_coords_file.pl  |\
    $MYPATH/reconcile_matches.pl $PREFIX.gap_coordinates.txt |\
    $MYPATH/output_reconciled_scaffolds.pl $HYB_CTG.broken|\
    ufasta format > $REF_CHR.$HYB_CTG.reconciled.fa && touch $PREFIX.place_extra.success
  else
    cat $PREFIX.reconciled.txt | $MYPATH/output_reconciled_scaffolds.pl $HYB_CTG.broken|\
    ufasta format > $REF_CHR.$HYB_CTG.reconciled.fa && touch $PREFIX.place_extra.success
  fi
fi

if [ -e $PREFIX.place_extra.success ];then
  log "Success! Final scaffold are in $REF_CHR.$HYB_CTG.reconciled.fa"
fi

#this line resets all gaps to 100
#perl -ane '{if($F[0] =~ /^>/){print;}else{@f=split(/N+/,$F[0]); print join("N"x100,@f),"\n"}}'  |\

