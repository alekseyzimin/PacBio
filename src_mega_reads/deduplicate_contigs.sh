#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
ASM_DIR=$1;
ASM_PREFIX=$2;
NUM_THREADS=$3;
PLOIDY=$4;

if [ $PLOIDY -gt 1 ];then
MERGE_LEN=40000
HAP_SIM_RATE=90
else
MERGE_LEN=10000
HAP_SIM_RATE=95
fi

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

set -e

#here we map the contigs against themselves to figure out which ones are redundant
perl -ane 'BEGIN{$seq=""}{if($F[0] =~ /^\>/){if(length($seq)>0){print length($seq)," $rn $seq\n";}$seq="";$rn=$F[0];}else{$seq.=$F[0];}}END{length($seq)," $rn $seq\n";}' $ASM_DIR/9-terminator/$ASM_PREFIX.scf.fasta | sort -S 50% -nrk1 | awk '{print $2"\n"$3}' >  $ASM_DIR/scaffolds.ref.fa && \
nucmer -t $NUM_THREADS --batch $(($(stat -c%s "$ASM_DIR/scaffolds.ref.fa")/25)) -l 31 -b 100 -p $ASM_DIR/sasm_to_sasm  $ASM_DIR/scaffolds.ref.fa  $ASM_DIR/9-terminator/$ASM_PREFIX.scf.fasta && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/sasm_to_sasm.delta > $ASM_DIR/sasm_to_sasm.noself.delta && \
parallel_delta-filter.sh $ASM_DIR/sasm_to_sasm.noself -q $NUM_THREADS && \
show-coords -lcHr -I $HAP_SIM_RATE $ASM_DIR/sasm_to_sasm.noself.fdelta | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl $MERGE_LEN | perl -ane '{$cov{$F[-1]}+=$F[15] if($F[15]>=10);}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>90);}}' > $ASM_DIR/sduplicates.txt && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/sasm_to_sasm.delta| show-coords -lcH /dev/stdin | awk '{if($12>$13 && $10>int("'$HAP_SIM_RATE'") && $16>90) print $NF}' >> $ASM_DIR/sduplicates.txt && \
ufasta extract -v -f $ASM_DIR/sduplicates.txt $ASM_DIR/9-terminator/$ASM_PREFIX.scf.fasta > $ASM_DIR/dedup.$ASM_PREFIX.scf.fasta
#echo "Final scaffold sequences are in $ASM_DIR/dedup.$ASM_PREFIX.scf.fasta"
