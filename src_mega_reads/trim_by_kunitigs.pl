#!/usr/bin/env perl
#this code trims mega-reads based on the unique k-unitigs on the ends of super-reads
my $superReadNamesFile=$ARGV[0];
my $kUnitigLengthsFile=$ARGV[1];
my $kmer=1000000;

open(FILE,$superReadNamesFile);
while($line=<FILE>){
chomp($line);
@f=split(/_/,$line);
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
  $trim_ku{$k}=1 if($sku{$k}==1&&not(defined($mku{$k}))&&not(defined($eku{$k})));
}
foreach $k(keys %eku){
  $trim_ku{$k}=1 if($eku{$k}==1&&not(defined($mku{$k}))&&not(defined($sku{$k})));
}

while($line=<STDIN>){
chomp($line);
@ff=split(/\s+/,$line);
@f=split(/_/,$ff[2]);
$start_trim=0;
$end_trim=0;
if($ff[1] eq "F"){
$start_trim=$len{substr($f[0],0,-1)}-$kmer if(defined($trim_ku{substr($f[0],0,-1)}));
$end_trim=$len{substr($f[-1],0,-1)}-$kmer if(defined($trim_ku{substr($f[-1],0,-1)}));
}else{
$end_trim=$len{substr($f[0],0,-1)}-$kmer if(defined($trim_ku{substr($f[0],0,-1)}));
$start_trim=$len{substr($f[-1],0,-1)}-$kmer if(defined($trim_ku{substr($f[-1],0,-1)}));
}
print "$ff[0] $start_trim $end_trim ",substr($f[0],0,-1)," ",substr($f[-1],0,-1),"\n";
}


