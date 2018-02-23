#!/usr/bin/env perl
#this code replaces consensus of one assembly with the other based on the nucmer matches -- consensus of the reference is replaced with the query
#
my $ref_contigs=$ARGV[0];
my $qry_contigs=$ARGV[1];

my $ctg="",$seq="";
my %rseq =(),%qseq=();
open(FILE,$ref_contigs);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $rseq{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $rseq{$ctg}=$seq;
}

$ctg="",$seq="";
open(FILE,$qry_contigs);
while($line=<FILE>){
  chomp($line);   
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $qseq{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $qseq{$ctg}=$seq;
}



#now read in the coords file
my $subseq;
my $last_offset=0;
my $offset=0;
while($line=<STDIN>){
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  next if($f[0]<=$last_offset);
  #$rseq{$f[-2]}="N"x$f[11] if(not(defined($rseq{$f[-2]})));
  next if(not(defined($rseq{$f[-2]})));
  die("query sequence $f[-1] not found") unless(defined($qseq{$f[-1]}));
  $subseq="";
  if($f[3]<$f[4]){
    $subseq=substr($qseq{$f[-1]},$f[3]-1,$f[4]-$f[3]+1);
  }else{
    $subseq=substr($qseq{$f[-1]},$f[4]-1,$f[3]-$f[4]+1);
    $subseq=reverse($subseq);
    $subseq=~tr/ACGTNacgtn/TGCAntgcan/;
  }
  $f[0]+=$offset;
  $f[1]+=$offset;
  $rseq{$f[-2]}=substr($rseq{$f[-2]},0,$f[0]-1).$subseq.substr($rseq{$f[-2]},$f[1]);
  $offset+=(length($subseq)-($f[1]-$f[0])-1);
  $last_offset=$f[1];
}

foreach $c(keys %rseq){
  print ">$c\n$rseq{$c}\n";
}
