#!/usr/bin/env perl
#
#this code finds a set of non-overlapping matches for a mega-read
#
#
#assume input sorted by mega-read size
$exedir=$ARGV[1];
open(FILE,$ARGV[0]);
$line=<FILE>;
($pbn,$pbseq)=split(/\s+/,$line);
close(FILE);
@interval_bgn=();
@interval_end=();
@mr=();
@nsr=();
@mrseq=();
$max_overlap=100;

while($line=<STDIN>){
chomp($line);
($nsr,$bgn,$end,$len,$mr,$mrseq)=split(/\s+/,$line);
$coords=`$exedir/swalign $pbseq $mrseq`;#refine coordinates with smith-waterman
@c=split(/\s+/,$coords);
$bgn=$c[2]-$c[0];
$end=$c[3]-$c[0];

for($i=0;$i<=$#interval_bgn;$i++){
if($bgn>=$interval_bgn[$i] && $bgn<=$interval_end[$i]){
$bgn_inside=1;
}else{
$bgn_inside=0;
}
if($end>$interval_bgn[$i] && $end<$interval_end[$i]){
$end_inside=1;
}else{
$end_inside=0;
}
last if($bgn_inside==1 && $end_inside==1);#fully contained
next if($bgn_inside==0 && $end_inside==0);#non-overlapping this interlal
if($bgn_inside==1){#check the overlaps
last if($interval_end[$i]-$bgn>$max_overlap);#overlap bigger -- contained
}else{
last if($end-$interval_bgn[$i]>$max_overlap);
}
}

if($i>$#interval_bgn){#not contained anywhere
push(@interval_bgn,$bgn);
push(@interval_end,$end);
push(@mr,$mr);
push(@nsr,$nsr);
push(@mrseq,$mrseq);
}
}

for($i=0;$i<=$#interval_bgn;$i++){
print "$nsr[$i] $interval_bgn[$i] $interval_end[$i] $mr[$i] $mrseq[$i]\n";
}

