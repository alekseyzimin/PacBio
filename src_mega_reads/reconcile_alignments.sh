#!/bin/bash
set -e
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

COORDS_FILE=$1;
SEQ_FILE=$2;

merge_matches_and_tile_coords_file.pl 1000000 < $COORDS_FILE | reconcile_matches.pl 
#| output_reconciled_scaffolds.pl $SEQ_FILE
