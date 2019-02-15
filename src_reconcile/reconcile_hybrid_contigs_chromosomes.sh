#!/bin/bash
#this code aims at reconciling the hybrid contigs and the chromosomes of the previously produces assembly
#arguments are: reference chromosomes, hybrid contigs, hybrid posmap (frgctg), filtered delta-file of alignments of ref to hyb
#MUST HAVE MaSURCA bin on the PATH
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
#low coverage threshold for breaking
COV_THRESH=3
REP_COV_THRESH=30

function error_exit {
    echo "$1" >&2
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
        -q|--query)
            QRY="$2"
            shift
            ;;
        -i|--identity)
            IDENTITY="$2"
            shift
            ;;
        -p|--posmap)
            POSMAP="$2"
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
        -r|--reference)
            REF="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage: reconcile_hybrid_contigs_chromosomes.sh -r <reference genome> -q <assembly to be scaffolded with the reference> -t <number of threads> -s <minimum sequence similarity percentage> -m <merge polishing sequence alignments slack (advanced)> -v <verbose> -p <posmap file for the assembly> -cl <coverage threshold for splitting at misassemblies, default 3> -ch <repeat coverage threshold for splitting at misassemblies, default 30>" 
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
HYB_CTG=`basename $QRY`
HYB_POS=$POSMAP

if [ ! -s $REF_CHR.$HYB_CTG.delta ];then
nucmer -t $NUM_THREADS -p $REF_CHR.$HYB_CTG -c 200 $REF $QRY
touch .rerun
fi

if [ ! -s $REF_CHR.$HYB_CTG.1.delta ] || [ -e .rerun ];then
delta-filter -1 -i $IDENTITY -o 20 $REF_CHR.$HYB_CTG.delta >$REF_CHR.$HYB_CTG.1.delta
touch .rerun
fi

#compute coverage from the posmap file
if [ ! -s $HYB_POS.coverage ];then 
awk '{print $1" "$2" "$3"\n"$1" "$2" "$4}' $HYB_POS | grep -v F |grep -v R | sort -nk2 -k3n -S 10% | compute_coverage.pl > $HYB_POS.coverage
touch .rerun
fi

if [ ! -s  gap_coordinates.txt ];then
splitFileAtNs $REF 1 > $REF_CHR.split.fasta
perl -ane '{$h{substr($F[1],3)}=$F[0]}END{while($line=<STDIN>){chomp($line);@f=split(/\s+/,$line);print "$f[0] $h{$f[1]} ",$f[2]+1," $f[3] $f[4]\n";}}' scaffNameTranslations.txt < genome.posmap.ctgscf | awk 'BEGIN{pg=0}{print $2" "pg" "$3;pg=$4}' > gap_coordinates.txt
touch .rerun
fi

if [ ! -s $REF_CHR.$HYB_CTG.1.coords ] || [ -e .rerun ];then
show-coords -lcHr $REF_CHR.$HYB_CTG.1.delta | \
merge_matches_and_tile_coords_file.pl $MERGE | \
merge_matches_and_tile_coords_file.pl $(($MERGE/10)) | \
awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
awk '{if($16>5 || $7>5000 ) print $0}' > $REF_CHR.$HYB_CTG.1.coords
touch .rerun
fi

#find split contigs
awk '{if($4<$5) print $4" "$5" "($4+$5)/2" "$NF" "$13; else print $5" "$4" "($4+$5)/2" "$NF" "$13;}' $REF_CHR.$HYB_CTG.1.coords| \
sort -k4 -k1n | \
uniq -D -f 3 | \
awk '{if($NF != prev){offset=$2;prev=$NF;print $0}else if($2>offset){print $0;offset=$2;}else{print "contained "$0}}' > $REF_CHR.$HYB_CTG.1.coords.split_contain

#create breaks files
grep -v ^contained $REF_CHR.$HYB_CTG.1.coords.split_contain | uniq -D -f 3 | awk '{if($4==ctg){if($1>5000 && $1<$NF-5000) print "alnbreak "substr($4,4)" "$1" 0"}else{ctg=$4;if($2>5000 && $2<$NF-5000) print "alnbreak "substr($4,4)" "$2" 0"}}' > $REF_CHR.$HYB_CTG.1.coords.breaks

if [ ! -s $HYB_POS.coverage.w_breaks ] || [ -e .rerun ];then
cat $REF_CHR.$HYB_CTG.1.coords.breaks $HYB_POS.coverage  | sort -nk2 -k3n -S 10% > $HYB_POS.coverage.w_breaks
touch .rerun
fi

grep -C 50 break $HYB_POS.coverage.w_breaks  | evaluate_splits.pl <(ufasta sizes -H $QRY | awk '{print substr($1,4)" "$2}') | sort -nk3 -S 10% >  $HYB_POS.coverage.w_breaks.validated

if [ ! -s $HYB_CTG.broken ] || [ -e .rerun ];then
break_contigs.pl <(grep -v "end" $HYB_POS.coverage.w_breaks.validated |awk '{if($4<=int("'$COV_THRESH'") || ($1="repeat" && $4>=int("'$REP_COV_THRESH'"))) print $0}') < $QRY > $HYB_CTG.broken
touch .rerun
fi

#now we re-align the broken contigs to the reference
if [ ! -s $REF_CHR.$HYB_CTG.broken.delta ] || [ -e .rerun ];then
nucmer -t $NUM_THREADS -p $REF_CHR.$HYB_CTG.broken -c 200  $REF $HYB_CTG.broken
touch .rerun
fi
if [ ! -s $REF_CHR.$HYB_CTG.broken.1.delta ] || [ -e .rerun ];then
delta-filter -1 -o 20 -i $IDENTITY $REF_CHR.$HYB_CTG.broken.delta > $REF_CHR.$HYB_CTG.broken.1.delta
touch .rerun
fi

rm -rf .rerun

#now we merge/rebuild chromosomes
show-coords -lcHr $REF_CHR.$HYB_CTG.broken.1.delta | \
merge_matches_and_tile_coords_file.pl $MERGE | \
merge_matches_and_tile_coords_file.pl $(($MERGE/10)) | \
awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=$2;last_scf=$18} if($2>last_end-10000) {print $0; last_end=$2}}' | \
awk '{if($16>5 || $7>5000 ) print $0}' > $REF_CHR.$HYB_CTG.broken.1.coords

# here we split everything so the contigs are "perfect"
cat $REF_CHR.$HYB_CTG.broken.1.coords  |  extract_single_best_match_coords_file.pl  |reconcile_matches.pl gap_coordinates.txt  > reconciled_coords.txt

cat reconciled_coords.txt  | output_reconciled_scaffolds.pl $HYB_CTG.broken| perl -ane '{if($F[0] =~ /^>/){print;}else{@f=split(/N+/,$F[0]); print join("N"x100,@f),"\n"}}'  | ufasta format > $REF_CHR.$HYB_CTG.reconciled.fa

