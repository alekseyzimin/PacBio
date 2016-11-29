#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
CA_PATH=$1
CA=$2;
ASM_PREFIX=$3;
NUM_THREADS=$4;
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$CA_PATH:$MYPATH:$PATH;

set -e

NUM_SUPER_READS=`cat mr.*.frg | grep -c --text '^{FRG' `
(cd $CA;
tigStore -g $ASM_PREFIX.gkpStore -t $ASM_PREFIX.tigStore 5 -d layout -U | \
tr -d '-' | \
awk 'BEGIN{print ">unique unitigs"}{if($1 == "cns"){seq=$2}else if($1 == "data.unitig_coverage_stat" && $2>=5){print seq"N"}}' | \
jellyfish count -L 2 -C -m 22 -s $ESTIMATED_GENOME_SIZE -t $NUM_THREADS -o unitig_mers /dev/fd/0;

cat <(overlapStore -b 1 -e $NUM_SUPER_READS -d genome.ovlStore  | awk '{if($1<$2 && $1<'$NUM_SUPER_READS' && $2<'$NUM_SUPER_READS') print $0}' | \
filter_overlap_file -t $NUM_THREADS <(gatekeeper  -dumpfragments -withsequence genome.gkpStore| grep -P '^fragmentIdent|^fragmentSequence' | \
perl -ane 'BEGIN{$flag=1}{if($flag){print ">";}print "$F[2]\n";$flag=1-$flag;}') unitig_mers /dev/fd/0) <(overlapStore -d genome.ovlStore | \
awk '{if($1<$2 && ($1>='$NUM_SUPER_READS' || $2>='$NUM_SUPER_READS')) print $1" "$2" "$3" "$4" "$5" "$6" "$7}')  |convertOverlap -ovl |gzip -c > overlaps.ovb.gz;
cd ..;
)
