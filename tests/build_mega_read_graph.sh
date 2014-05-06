#!/bin/bash

FILENAME=$1  #filename of a megareads file for matches of SINGLE pacbio read to the super-reads
NAMESEQFILE=$2 #file containing pb read "name sequence"
KMER=$3
let KM1=$KMER-1
EXEPATH=`dirname $0`

sort -nrk3,3 ${FILENAME} > $FILENAME.megareads

#wc -l $FILENAME.megareads
if [ -s $FILENAME.megareads ];then
/home/alekseyz/myprogs/MaSuRCA/build/inst/bin/createFastaSuperReadSequences work1 <(awk '{print "1 "$6}' $FILENAME.megareads) -seqdiffmax 0 -min-ovl-len $KM1 -minreadsinsuperread 1  -good-sr-filename $FILENAME.megareads.names  -kunitigsfile /genome7/raid/alekseyz/PB_ScerW303/assembly_k${KMER}/guillaumeKUnitigsAtLeast32bases_all.fasta -good-sequence-output-file $FILENAME.megareads.fa -super-read-name-and-lengths-file $FILENAME.megareads.sizes -rename-super-reads  2>/dev/null

awk -F ',' '{print ">"$1"\n"$2}'  $NAMESEQFILE > $NAMESEQFILE.fa
nucmer  -d 0.3 -f -g 300 -l 15 -b 1000 -p $FILENAME $NAMESEQFILE.fa $FILENAME.megareads.fa 1>/dev/null 2>&1
delta-filter -1 $FILENAME.delta > $FILENAME.f.delta
#show-coords -lcHr -I 75 $FILENAME.delta | $EXEPATH/extract_best_match_coords.pl > $FILENAME.f.ncoords
show-coords -lcHr -I 75 $FILENAME.f.delta | /home/alekseyz/myprogs/merge_matches_coords_file.pl > $FILENAME.f.ncoords

perl -e '{
open(FILE,$ARGV[0]);
$n=0;
while($line=<FILE>){
	chomp($line);
	($srname,$srl)=split(/\s+/,$line);
	$srn{$srname}=$n;
	push(@srnames,$srname);
	$n++;
}
open(FILE,$ARGV[1]);
while($line=<FILE>){
	chomp($line);
	$rn=substr($line,1);
	$line=<FILE>;
	chomp($line);
	push(@seq,$line);
}
open(FILE,$ARGV[2]);
while($line=<FILE>){
	chomp($line);
	@f=split(/\s+/,$line);
	$scores[$srn{$f[-1]}]=$f[2];
	$coords_pb[$srn{$f[-1]}]="$f[0] $f[1]";
	$pbn[$srn{$f[-1]}]=$f[4];
}
while($line=<STDIN>){
	chomp($line);
	$line=~s/^\s+//;
	@f=split(/\s+/,$line);
	@c=split(/\s+/,$coords_pb[$f[-1]]);
        next if(($f[0]>$c[1] && $f[1]>$c[1])||($f[0]<$c[0] && $f[1]<$c[0]));
	$score=int($f[7]*$f[9]/100);
	if($score>$scores[$f[-1]]){
	$scores[$f[-1]]=$score;
	$outline[$f[-1]]="$f[0] $f[1] $f[3] $f[4] ".($f[1]-$f[0])." $pbn[$f[-1]] $seq[$f[-1]] $srnames[$f[-1]]\n";
	}

}
foreach $l(@outline){
print $l;
}
}' $FILENAME.megareads.sizes $FILENAME.megareads.fa $FILENAME.megareads < $FILENAME.f.ncoords  | sort -nrk5,5 > $FILENAME.megareads.wseq


$EXEPATH/reconciliate_mega_reads.pl < $FILENAME.megareads.wseq | sort -nk1,1 > $FILENAME.megareads.sorted.wseq
fi
