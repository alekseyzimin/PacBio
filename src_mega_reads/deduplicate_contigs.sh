#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

ASM_FASTA=$1;

nucmer -l 31 -c 400 -b 400 -p asm_to_asm $ASM_FASTA $ASM_FASTA
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' asm_to_asm.delta| delta-filter -q -o 20 /dev/stdin|show-coords -lcHr /dev/stdin | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file.pl 100000 | awk '{if($10>95 && $16>90) print $NF;}' > $ASM_FASTA.duplicates.txt
ufasta extract -f $ASM_FASTA.duplicates.txt $ASM_FASTA > $ASM_FASTA.alt
ufasta extract -v -f $ASM_FASTA.duplicates.txt $ASM_FASTA > $ASM_FASTA.main
