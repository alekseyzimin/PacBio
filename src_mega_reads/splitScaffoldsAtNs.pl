#!/usr/bin/env perl
#this script splits scaffolds at Ns
while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
  if(length($seq)>0){
      $seq=~s/n/N/g;
      @f=split(/(N{1,})/,$seq); 
      my $n=1;
      foreach $c(@f){
        if(not($c=~/^N/) && length($c)>0){
          $start=$n;
          $end=$n+length($c)-1;
          print ">$rn.$end\n$c\n";
        }
        $n+=length($c);
      }
  }
  ($rn,$junk)=split(/\s+/,substr($line,1));
  $seq="";
  }else{
  $seq.=$line;
}
}
  $seq=~s/n/N/g;
  @f=split(/(N{1,})/,$seq);
  my $n=1;
  foreach $c(@f){
    if(not($c=~/^N/) && length($c)>0){
      $start=$n;
      $end=$n+length($c)-1;
      print ">$rn.$end\n$c\n";
    }
    $n+=length($c);
  }

