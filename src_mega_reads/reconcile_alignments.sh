#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

COORDS_FILE=$1;
SEQ_FILE=$2;

merge_matches_and_tile_coords_file.pl 1000000 < $COORDS_FILE | awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=0;last_scf=$18} if($2>last_end) print $0; last_end=$2}' |  reconcile_matches.pl 
#| output_reconciled_scaffolds.pl $SEQ_FILE
