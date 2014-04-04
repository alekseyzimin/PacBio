#!/usr/bin/env perl
#
#this code finds a set of non-overlapping matches for a mega-read
#
#
#assume input sorted by mega-read size

@interval_bgn=();
@interval_end=();
@mr=();
@nsr=();
@mrseq=();
$max_overlap=69;

while($line=<STDIN>){
chomp($line);
($pbgn,$pend,$mbgn,$mend,$qlt,$pb,$mrseq)=split(/\s+/,$line);
$mrseq=substr($mrseq,$mbgn-1,$mend-$mbgn+1);
$bgn=$pbgn;
$end=$pend;
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
push(@mrseq,$mrseq);
}
}

for($i=0;$i<=$#interval_bgn;$i++){
print "$interval_bgn[$i] $interval_end[$i] $pb $mrseq[$i]\n";
}

