#!/usr/bin/env perl
#
#this code finds a set of non-overlapping matches for a mega-read
#
#
#assume input sorted by mega-read size

@interval_bgn=();
@interval_end=();
@interval_g_bgn=();
@interval_g_end=();

@mr=();
@nsr=();
@mrseq=();
@mrname=();
#max overlap is in percentage of size
$max_overlap_pct=$ARGV[0];

while($line=<STDIN>){
chomp($line);
($pbgn,$pend,$mbgn,$mend,$qlt,$pb,$mrseq,$mrname)=split(/\s+/,$line);
$mrseq=substr($mrseq,$mbgn-1,$mend-$mbgn+1);
$max_overlap=$max_overlap_pct*length($mrseq)/100;
$bgn=$pbgn;
$end=$pend;
$overlap=0;
for($i=0;$i<=$#interval_g_bgn;$i++){

last if($bgn >= $interval_g_bgn[$i] && $end <= $interval_g_end[$i]);#contained
last if($bgn < $interval_g_bgn[$i] && $end > $interval_g_end[$i]);#containing -- happens on rare occasions

$bgn_inside=0;
$end_inside=0;
$bgn_inside=1 if($bgn>=$interval_g_bgn[$i] && $bgn<=$interval_g_end[$i]);
$end_inside=1 if($end>=$interval_g_bgn[$i] && $end<=$interval_g_end[$i]);
next if($bgn_inside==0 && $end_inside==0);#non-overlapping this interlal

#we get here if we have an overlap
if($bgn_inside==1){#check the overlaps
last if($interval_g_end[$i]-$bgn>$max_overlap);#overlap bigger -- contained
$interval_g_end[$i]=$end;
$overlap=1;
}else{
last if($end-$interval_g_bgn[$i]>$max_overlap);
$interval_g_bgn[$i]=$bgn;
$overlap=1;
}
}
if($i>$#interval_g_bgn){#not contained anywhere
if($overlap==0){
push(@interval_g_bgn,$bgn);
push(@interval_g_end,$end);
}
push(@interval_bgn,$bgn);
push(@interval_end,$end);
push(@mrseq,$mrseq);
push(@mrnames,$mrname);
}
}

for($i=0;$i<=$#interval_bgn;$i++){
print "$interval_bgn[$i] $interval_end[$i] $pb $mrseq[$i] $mrnames[$i]\n";
}

