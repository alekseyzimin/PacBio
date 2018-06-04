#!/usr/bin/env perl
#this code trims mega-reads given the trim points
my $trimsFile=$ARGV[0];
open(FILE,$trimsFile);
while($line=<FILE>){
  @f=split(/\s+/,$line);
  $frontTrim{$f[0]}=$f[1];
  $backTrim{$f[0]}=$f[2];
}

while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
    @f=split(/\s+/,$line);
    $readName=substr($f[0],1);
  }else{
    if($backTrim{$readName}>0){
      $seq=substr($line,$frontTrim{$readName},-$backTrim{$readName}),"\n";
    }else{
      $seq=substr($line,$frontTrim{$readName}),"\n";
    }
    print ">$readName\n$seq\n" if(length($seq)>=500);
  }
}
