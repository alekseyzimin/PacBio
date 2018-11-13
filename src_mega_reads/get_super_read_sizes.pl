#!/usr/bin/env perl
#this code creates super-reads size file from readPositionsInSuperReads
my $kuSizeFile=$ARGV[0];
my $line,$minkusize=10000000;
open(FILE,$kuSizeFile); #output of ufasta sizes -H on kunitigs file
while($line=<FILE>){
  chomp($line);
  my ($ku,$sz)=split(/\s+/,$line);
  $kusize[$ku]=$sz;
  $minkusize=$sz if($sz< $minkusize);
  }
$minkusize--;#need to know k-1
#now we read in the super-read names WITHOUT F and R  and compute the super-read sizes
for(my $i=0;$i<=$#kusize;$i++){
$kusize[$i]-=$minkusize;
}
while($line=<STDIN>){
  chomp($line);
  @srname=split(/_/,$line);
  my $srsize=$kusize[substr($srname[0],0,-1)]+$minkusize;
  
  for(my $i=1;$i<=$#srname; $i++){
    $srsize+=$kusize[substr($srname[$i],0,-1)];
    }
    print "$line $srsize\n";
}

