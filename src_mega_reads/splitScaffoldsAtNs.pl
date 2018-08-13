#!/usr/bin/env perl

while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
  if(length($seq)>0){
      @f=split(/(N{1,})/,uc($seq)); 
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
  $rn=substr($line,1);
  $seq="";
  }else{
  $seq.=$line;
}
}

  @f=split(/(N{1,})/,uc($seq));
  my $n=1;
  foreach $c(@f){
    if(not($c=~/^N/) && length($c)>0){
      $start=$n;
      $end=$n+length($c)-1;
      print ">$rn.$end\n$c\n";
    }
    $n+=length($c);
  }

