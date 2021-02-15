#!/bin/sh
FILE=$1;
NUM=$2;
perl -ane '{
if($F[0]=~/^>/){
  if(length($seq)>0){
      $seq=~s/n/N/g;
      @f=split(/(N{'$NUM',})/,$seq); 
      my $n=1;
      foreach $c(@f){
        if(not($c=~/^N/) && length($c)>0){
          $start=$n;
          $end=$n+length($c)-1;
          print ">$rn:$start-$end\n$c\n";
        }
        $n+=length($c);
      }
  }
  $rn=substr($F[0],1);
  $seq="";
  }else{
  $seq.=$F[0];
}
}END{
  $seq=~s/n/N/g;
  @f=split(/(N{'$NUM',})/,$seq);
  my $n=1;
  foreach $c(@f){
    if(not($c=~/^N/) && length($c)>0){
      $start=$n;
      $end=$n+length($c)-1;
      print ">$rn:$start-$end\n$c\n";
    }
    $n+=length($c);
  }
}'  $FILE
