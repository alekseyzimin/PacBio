#!/bin/bash

FILENAME=$1  #filename of a megareads file for matches of SINGLE pacbio read to the super-reads
NAMESEQFILE=$2 #file containing pb read "name sequence"
KMER=$3
let KM1=$KMER-1
EXEPATH=`dirname $0`

sort -nrk3,3 ${FILENAME} > $FILENAME.megareads
rm -f $FILENAME.megareads.sorted.wseq
touch $FILENAME.megareads.sorted.wseq

#wc -l $FILENAME.megareads
if [ -s $FILENAME.megareads ];then
/home/alekseyz/myprogs/masurca-devel/build/inst/bin/createFastaSuperReadSequences work1 <(awk '{print "1 "$6}' $FILENAME.megareads| head -n 25) -seqdiffmax 0 -min-ovl-len $KM1 -minreadsinsuperread 1  -good-sr-filename $FILENAME.megareads.names  -kunitigsfile /genome7/raid/alekseyz/PB_ScerW303/assembly_k${KMER}/guillaumeKUnitigsAtLeast32bases_all.fasta -good-sequence-output-file $FILENAME.megareads.all.fa -super-read-name-and-lengths-file $FILENAME.megareads.sizes -rename-super-reads  2>/dev/null

perl -e 'while($line=<STDIN>){$s=$line;$line=<STDIN>;print "$s$line" if(length($line)>125);}' < $FILENAME.megareads.all.fa > $FILENAME.megareads.fa
if [ ! -s $FILENAME.megareads.fa ];then exit; fi
awk -F ',' '{print ">"$1"\n"$2}'  $NAMESEQFILE > $NAMESEQFILE.fa
nucmer  -d 0.2 -f -g 200 -l 15 -b 200 -p $FILENAME $NAMESEQFILE.fa $FILENAME.megareads.fa 1>/dev/null 2>&1
if [ ! -s $FILENAME.delta ];then exit; fi
delta-filter -g -o 20 $FILENAME.delta > $FILENAME.f.delta
show-coords -lcHr $FILENAME.f.delta  |sort -nk19 -k1n |~/myprogs/merge_matches_coords_file.pl 1000 > $FILENAME.f.ncoords
if [ ! -s $FILENAME.f.ncoords ];then exit; fi

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
	next if (not(defined($srn{$f[5]})));
	$scores[$srn{$f[5]}]=$f[2];
	$coords_pb[$srn{$f[5]}]="$f[0] $f[1]";
	$pbn[$srn{$f[5]}]=$f[4];
}
while($line=<STDIN>){
	chomp($line);
	$line=~s/^\s+//;
	@f=split(/\s+/,$line);
	@c=split(/\s+/,$coords_pb[$f[-1]]);
        next if(($f[0]>$c[1] && $f[1]>$c[1])||($f[0]<$c[0] && $f[1]<$c[0]));
	$outline[$f[-1]].="$f[0] $f[1] $f[3] $f[4] ".$scores[$f[-1]]." $pbn[$f[-1]] $seq[$f[-1]] $srnames[$f[-1]]\n";
}
foreach $l(@outline){
print $l;
}
}' $FILENAME.megareads.sizes $FILENAME.megareads.fa $FILENAME.megareads < $FILENAME.f.ncoords  | sort -nrk5,5 > $FILENAME.megareads.wseq


$EXEPATH/reconciliate_mega_reads.pl 20 < $FILENAME.megareads.wseq | sort -nk1,1 > $FILENAME.megareads.sorted.wseq
fi
