#!/usr/bin/env perl
#this code trims mega-reads based on the unique k-unitigs on the ends of super-reads
my $readPlacementFile=$ARGV[0];
my $megaReadsFile=$ARGV[1];
my $superReadNamesSizesFile=$ARGV[2];
my $kUnitigLengthsFile=$ARGV[3];
my $kmer=1000000;
my @mr_sizes;
my @mr_names;
my %groups;

open(FILE,$megaReadsFile);
while($line=<FILE>){
  if(substr($line,0,1) eq ">"){
    chomp($line);
    push(@mr_names,substr($line,1));
  }else{
    push(@mr_sizes,length($line));
  }
}

my $numKUnitigs=0;
open(FILE,$kUnitigLengthsFile);
while($line=<FILE>){
  chomp($line);
  my ($kunitig,$kulen)=split(/\s+/,$line);
  $len[$kunitig]=$kulen;
  $kmer=$kulen if($kmer>$kulen);
  $numKUnitigs++;
}

my @sku=(0)x$numKUnitigs;
my @mku=(0)x$numKUnitigs;
my @eku=(0)x$numKUnitigs;
my @trim_ku=(0)x$numKUnitigs;

open(FILE,$superReadNamesSizesFile);
while($line=<FILE>){
  chomp($line);
  my ($name,$length)=split(/\s+/,$line);
  @f=split(/_/,$name);
  next if($#f<2);
  $sku[substr($f[0],0,-1)]++;
  for($i=1;$i<$#f;$i++){
    $mku[substr($f[$i],0,-1)]++;
  }
  $eku[substr($f[-1],0,-1)]++;
}
close(FILE);

$kmer--;

for($k=0;$k<$numKUnitigs;$k++){
  $trim_ku[$k]=1 if($sku[$k]==1 && $mku[$k]==0 && $eku[$k]==0); 
  $trim_ku[$k]=1 if($eku[$k]==1 && $mku[$k]==0 && $sku[$k]==0);
}

open(FILE,$readPlacementFile);
while($line=<FILE>){
  chomp($line);
  my ($read,$sread,$pos,$ori,$code)=split(/\s+/,$line);
  my @f=split(/_/,$sread);
  my $start_trim=0;
  my $end_trim=0;
  my $firstKU=substr($f[0],0,-1);
  my $lastKU=substr($f[-1],0,-1);
  if($ori eq "F"){
    $start_trim=$len[$firstKU]-$kmer if($trim_ku[$firstKU]);
    $end_trim=$len[$lastKU]-$kmer if($trim_ku[$lastKU]);
  }else{
    $end_trim=$len[$firstKU]-$kmer if($trim_ku[$firstKU]);
    $start_trim=$len[$lastKU]-$kmer if($trim_ku[$lastKU]);
  }
  print $mr_names[int(substr($read,2)/2)]," $start_trim $end_trim $firstKU $lastKU\n";
}


