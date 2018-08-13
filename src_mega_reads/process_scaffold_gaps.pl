#!/usr/bin/env perl
#this code extracts mate pairs that correctly span scaffold gaps
my $ctgscf=$ARGV[0];
my $frgctg=$ARGV[1];
my $DEBUG=0;
#reading in contig pairs
open(FILE,$ctgscf);
my $prevscf="";
my $prevctg="";
while($line=<FILE>){
  chomp($line);
  my @f=split(/\s+/,$line);
  if($f[1] eq $prevscf){
    $ctp{"$prevctg$prevdir $f[0]$f[4]"}="$f[1].$f[3]";
    $prevdir=~tr/fr/rf/;
    $f[4]=~tr/fr/rf/;
    $ctp{"$f[0]$f[4] $prevctg$prevdir"}="$f[1].$f[3]";
    $f[4]=~tr/fr/rf/;
  }
  $prevctg=$f[0];
  $prevdir=$f[4];
  $prevscf=$f[1];
}

if($DEBUG){
  foreach $k(keys %ctp){
    print "$k $ctp{$k}\n";
  }
}

open(FILE,$frgctg);
while($line=<FILE>){
  chomp($line);
  my @f=split(/\s+/,$line);
  my $mname=substr($f[0],0,length($f[0])-1);
  my $mdir=substr($f[0],length($f[0])-1);
  if($mdir eq "F"){
    if($f[4] eq "f"){
      $mctp{$mname}="$f[1]f".$mctp{$mname};
    }else{
      $mctp{$mname}="$f[1]r".$mctp{$mname};
    }
  }elsif($mdir eq "R"){
    if($f[4] eq "f"){
      $mctp{$mname}=$mctp{$mname}." $f[1]r";
    }else{
      $mctp{$mname}=$mctp{$mname}." $f[1]f";
    }
  }
}

foreach $k(keys %mctp){
  my ($rn,$ri)=split(/\./,$k);
  print "$rn scf$ctp{$mctp{$k}}\n" if(defined($ctp{$mctp{$k}}));
}

