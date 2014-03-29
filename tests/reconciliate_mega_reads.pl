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

while($line=<STDIN>){
chomp($line);
($nsr,$bgn,$end,$len,$mr)=split(/\s+/,$line);
for($i=0;$i<=$#interval_bgn;$i++){
last if($bgn>$interval_bgn[$i] && $bgn<$interval_end[$i]);
last if($end>$interval_bgn[$i] && $end<$interval_end[$i]);
}
if($i>$#interval_bgn){#not contained anywhere
push(@interval_bgn,$bgn);
push(@interval_end,$end);
push(@mr,$mr);
push(@nsr,$nsr);
}
}

for($i=0;$i<=$#interval_bgn;$i++){
print "$nsr[$i] $interval_bgn[$i] $interval_end[$i] $mr[$i]\n";
}

