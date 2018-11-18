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

open(FILE,$superReadNamesSizesFile);
while($line=<FILE>){
  chomp($line);
  my ($name,$length)=split(/\s+/,$line);
  @f=split(/_/,$name);
  next if($#f<2);
  $sku{substr($f[0],0,-1)}++;
  for($i=1;$i<$#f;$i++){
    $mku{substr($f[$i],0,-1)}++;
  }
  $eku{substr($f[-1],0,-1)}++;
}
close(FILE);

open(FILE,$kUnitigLengthsFile);
while($line=<FILE>){
  chomp($line);
  @ff=split(/\s+/,$line);
  $len{$ff[0]}=$ff[1];
  $kmer=$ff[1] if($kmer>$ff[1]);
}

$kmer--;

foreach $k(keys %sku){
  $trim_ku{$k}=1 if($sku{$k}==1 && not(defined($mku{$k}))&&not(defined($eku{$k})));
}
foreach $k(keys %eku){
  $trim_ku{$k}=1 if($eku{$k}==1 && not(defined($mku{$k}))&&not(defined($sku{$k})));
}

open(FILE,$readPlacementFile);
while($line=<FILE>){
  chomp($line);
  my ($read,$sread,$pos,$ori,$code)=split(/\s+/,$line);
  my @ff=($mr_names[int(substr($read,2)/2)],$ori,$sread);
  my @f=split(/_/,$ff[2]);
  my $start_trim=0;
  my $end_trim=0;
  if($ff[1] eq "F"){
    $start_trim=$len{substr($f[0],0,-1)}-$kmer if(defined($trim_ku{substr($f[0],0,-1)}));
    $end_trim=$len{substr($f[-1],0,-1)}-$kmer if(defined($trim_ku{substr($f[-1],0,-1)}));
  }else{
    $end_trim=$len{substr($f[0],0,-1)}-$kmer if(defined($trim_ku{substr($f[0],0,-1)}));
    $start_trim=$len{substr($f[-1],0,-1)}-$kmer if(defined($trim_ku{substr($f[-1],0,-1)}));
  }
  print "$ff[0] $start_trim $end_trim ",substr($f[0],0,-1)," ",substr($f[-1],0,-1),"\n";
}


