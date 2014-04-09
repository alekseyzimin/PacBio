#!/bin/bash

FILENAME=$1  #filename of a coords file for matches of SINGLE pacbio read to the super-reads
NAMESEQFILE=$2 #file containing pb read "name sequence"
EXEPATH=`dirname $0`

#will fail if there is a gap in pacBio read coverage

perl -ane '{$ms=$F[0]-$F[2];$me=$F[1]+$F[10]-$F[3]; print "$ms $me ",join(" ",@F[2..$#F]),"\n"}' $FILENAME | $EXEPATH/longest_path_overlap_graph -m 1 -H -d 0.05 | sort -nrk3,3 > $FILENAME.megareads
if [ -L /genome7/raid/alekseyz/PB_ScerW303/assembly/work1 ];then
ln -s /genome7/raid/alekseyz/PB_ScerW303/assembly/work1
fi
wc -l $FILENAME.megareads
if [ -s $FILENAME.megareads ];then
/home/alekseyz/myprogs/MaSuRCA/build/inst/bin/createFastaSuperReadSequences work1 <(awk '{print "1 "$6}' $FILENAME.megareads) -seqdiffmax 0 -min-ovl-len 69 -minreadsinsuperread 1  -good-sr-filename $FILENAME.megareads.names  -kunitigsfile /genome7/raid/alekseyz/PB_ScerW303/assembly/guillaumeKUnitigsAtLeast32bases_all.fasta -good-sequence-output-file $FILENAME.megareads.fa -super-read-name-and-lengths-file $FILENAME.megareads.sizes -rename-super-reads  2> work1/createFastaSuperReadSequences.errors.txt

awk -F ',' '{print ">"$1"\n"$2}'  $NAMESEQFILE > $NAMESEQFILE.fa
nucmer --maxmatch -d 0.25 -f -g 200 -l 13 -b 1000 -p $FILENAME $NAMESEQFILE.fa $FILENAME.megareads.fa 1>/dev/null 2>&1

#show-coords -lcHr -I 75 $FILENAME.delta | $EXEPATH/extract_best_match_coords.pl > $FILENAME.f.ncoords
show-coords -lcHr -I 75 $FILENAME.delta > $FILENAME.f.ncoords

perl -e '{
open(FILE,$ARGV[0]);
$n=0;
while($line=<FILE>){
	chomp($line);
	($srn,$srl)=split(/\s+/,$line);
	$srn{$srn}=$n;
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
	$outline[$f[-1]]="$f[0] $f[1] $f[3] $f[4] ".($f[1]-$f[0])." $pbn[$f[-1]] $seq[$f[-1]]\n";
	}

}
foreach $l(@outline){
print $l;
}
}' $FILENAME.megareads.sizes $FILENAME.megareads.fa $FILENAME.megareads < $FILENAME.f.ncoords  | sort -nrk5,5 > $FILENAME.megareads.wseq


$EXEPATH/reconciliate_mega_reads.pl < $FILENAME.megareads.wseq | sort -nk1,1 > $FILENAME.megareads.sorted.wseq
fi
