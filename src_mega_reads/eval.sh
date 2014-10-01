#!/bin/bash
COORDS=$1;
./create_fasta_frg.sh 200 $COORDS 70 pb10x.fasta;
cp $COORDS.200.fa $COORDS.200.n.fa;
nucmer -p $COORDS -l 51 ../finished_sequence/polished_assembly.fasta $COORDS.200.n.fa;
delta-filter -q -o 20 $COORDS.delta > $COORDS.q.delta;
show-coords -lcHr $COORDS.q.delta > $COORDS.q.coords;
awk '{if($1>500 && $2< $12-500){if($10>98&&$16>99){n+=$8}else {print $0" blasr"}m+=$8}}END{print n/m}' $COORDS.q.coords |wc -l
