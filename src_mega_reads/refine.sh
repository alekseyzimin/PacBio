#!/bin/bash
PACBIO=$1
COORDS=$2
INPUT=$3
KMER=$4
#set -o pipefail
#set -e

#exec > $INPUT.log 2>&1

refine_alignments.pl $PACBIO $INPUT < $INPUT | \
     delta-filter -r1 -o 50 /dev/stdin | \
     show-coords -lcHr /dev/stdin | \
     awk '{if($4<$5){print $18"/0_"$12" "$19" 0 0 0 "$10" "$4" "$5" "$13" "$1" "$2" "$12" 0"}else{print $18"/0_"$12" "$19+1" 0 0 0 "$10" "$13-$4+1" "$13-$5+1" "$13" "$1" "$2" "$12" 0"}}' > t.$INPUT.blasr.out && reconciliate_mega_reads.maximal.nucmer.pl 20 $KMER t.$INPUT.maximal_mr.fa t.$INPUT.maximal_mr.names <t.$INPUT.blasr.out  1>$COORDS$INPUT.all.txt.tmp 2>/dev/null && rm t.$INPUT.blasr.out t.$INPUT.maximal_mr.fa t.$INPUT.maximal_mr.names

