#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

COORDS_FILE=$1;
SEQ_FILE=$2;
REF_FILE=$3;
splitFileAtNs $REF_FILE 1 > /dev/null
perl -ane '{$h{substr($F[1],3)}=$F[0]}END{while($line=<STDIN>){chomp($line);@f=split(/\s+/,$line);print "$f[0] $h{$f[1]} ",$f[2]+1," $f[3] $f[4]\n";}}' scaffNameTranslations.txt < genome.posmap.ctgscf | awk 'BEGIN{pg=0}{print $2" "pg" "$3;pg=$4}' > gap_coordinates.txt
merge_matches_and_tile_coords_file.pl 500000 < $COORDS_FILE |merge_matches_and_tile_coords_file.pl 100000 | awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=0;last_scf=$18} if($2>last_end){ print $0; last_end=$2}}' |  awk '{if($16>5 || $7>5000) print $0}' > $COORDS_FILE.merged
cat $COORDS_FILE.merged  | reconcile_matches.pl gap_coordinates.txt | awk '{print $2" "$3" "$4}' | sort |uniq -d | awk '{print $1}' > duplicate_contigs.txt
cat $COORDS_FILE.merged  | reconcile_matches.pl gap_coordinates.txt duplicate_contigs.txt > reconciled_coords.txt
#| output_reconciled_scaffolds.pl $SEQ_FILE
