#!/usr/bin/env perl
my $pbseqfile =$ARGV[0];

open(FILE,$pbseqfile);
while($line=<FILE>){
  chomp($line);
  if($line =~ /^>/){
    @f=split(/\s+/,$line);
    $rn=substr($f[0],1);
  }else{
    $pbseq{$rn}.=$line;
  }
}

while($line=<STDIN>){
  chomp($line);
    if(substr($line,0,1) eq ">"){
    $rn=substr($line,1);
    die("read sequence for $rn not found") if(not(defined($pbseq{$rn})));
    print "$line $pbseq{$rn}\n";
    }else{
    print $line,"\n";
    }
}

