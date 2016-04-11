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

merge_matches_and_tile_coords_file.pl 500000 < $COORDS_FILE | awk '{if($16>1) print $0}' | merge_matches_and_tile_coords_file.pl 250000 | awk 'BEGIN{last_end=0;last_scf="";}{if($18 != last_scf){last_end=0;last_scf=$18} if($2>last_end){ print $0; last_end=$2}}' |  awk '{if($16>5 || $8>20000) print $0}' > $COORDS_FILE.merged

awk '{if($4<$5) print $4" "$5" "($4+$5)/2" "$NF; else print $5" "$4" "($4+$5)/2" "$NF;}'  $COORDS_FILE.merged |sort -k4 -k1n -S 10% |uniq -D -f 3 | awk '{if($NF != prev){offset=$2;prev=$NF;print $0}else if($2>offset){print $0;offset=$2;}else{print "contained "$0}}' > contig_split_contain.txt

grep contain contig_split_contain.txt | perl -ane '{$h{"$F[1] $F[2] $F[4]"}=1;}END{open(FILE,"'$COORDS_FILE'.merged");while($line=<FILE>){chomp($line);@f=split(/\s+/,$line);if($f[3]<$f[4]){$s=$f[3];$e=$f[4]}else{$s=$f[4];$e=$f[3]}; print $line,"\n" unless defined($h{"$s $e $f[-1]"});}}' | merge_matches_and_tile_coords_file.pl 250000 > $COORDS_FILE.nocontain

grep -v contain contig_split_contain.txt | uniq -D -f 3 | perl -ane 'BEGIN{$ctg="";}{if(not($F[3] eq $ctg)){if(@ctg){@f1=split(/\s+/,$ctg[0]);@f2=split(/\s+/,$ctg[1]);print "$f1[0] $f1[1] $f1[0] ",($f2[0]-$f1[1])/2," $f1[3]\n";for($i=1;$i<$#ctg;$i++){@f1=split(/\s+/,$ctg[$i-1]);@f2=split(/\s+/,$ctg[$i]);@f3=split(/\s+/,$ctg[$i+1]);print "$f2[0] $f2[1] ",($f2[0]-$f1[1])/2," ",($f3[0]-$f2[1])/2," $f2[3]\n";}@f1=split(/\s+/,$ctg[-2]);@f2=split(/\s+/,$ctg[-1]);print "$f2[0] $f2[1] ",($f2[0]-$f1[1])/2," 100000000 $f2[3]\n";}$ctg=$F[3];@ctg=();}push(@ctg,join(" ",@F))}END{if(@ctg){@f1=split(/\s+/,$ctg[0]);@f2=split(/\s+/,$ctg[1]);print "$f1[0] $f1[1] $f1[0] ",($f2[0]-$f1[1])/2," $f1[3]\n";for($i=1;$i<$#ctg;$i++){@f1=split(/\s+/,$ctg[$i-1]);@f2=split(/\s+/,$ctg[$i]);@f3=split(/\s+/,$ctg[$i+1]);print "$f2[0] $f2[1] ",($f2[0]-$f1[1])/2," ",($f3[0]-$f2[1])/2," $f2[3]\n";}@f1=split(/\s+/,$ctg[-2]);@f2=split(/\s+/,$ctg[-1]);print "$f2[0] $f2[1] ",($f2[0]-$f1[1])/2," 100000000 $f2[3]\n";}}' >restrict_extensions.txt

cat $COORDS_FILE.nocontain | reconcile_matches.pl gap_coordinates.txt restrict_extensions.txt > reconciled_coords.txt

cat reconciled_coords.txt | output_reconciled_scaffolds.pl $SEQ_FILE > $REF_FILE.reconciled.fa
