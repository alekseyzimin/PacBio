#!/usr/bin/env perl
use jellyfish;
use strict;
use warnings;

#ASSUMES one line per contig

my $noise_thresh=1;#kmers with counts above this threshold will be used
my $qfP1 = jellyfish::QueryMerFile->new($ARGV[0]);
my $qfP2 = jellyfish::QueryMerFile->new($ARGV[1]);
my $k = jellyfish::MerDNA::k();
my $mer = jellyfish::MerDNA->new;
my %countp1=();
my %countp2=();
my %countboth=();
my %count=();
my $ctg;
print STDERR "resolving haplotypes with k=$k\n";

while(my $line=<STDIN>){
  chomp($line);
  if($line =~ /^>/){
    my @f=split(/\s+/,$line);
    $ctg=substr($f[0],1);
    $countp1{$ctg}=0;
    $countp2{$ctg}=0;
    $countboth{$ctg}=0;
    $count{$ctg}=0;
  }else{
    for(my $i=0; $i<length($line)-$k+1; $i++){
      my $m=substr($line,$i,$k);
      $mer->set($m);
      $mer->canonicalize;
      my $count1=$qfP1->get($mer);
      my $count2=$qfP2->get($mer);
      $countp1{$ctg}++ if($count1>$noise_thresh);
      $countp2{$ctg}++ if($count2>$noise_thresh);
      $countboth{$ctg}++ if($count1>$noise_thresh && $count2>$noise_thresh);
      $count{$ctg}++;
    }
    print "$ctg ",$countp1{$ctg}-$countboth{$ctg}," ",$countp2{$ctg}-$countboth{$ctg}," $countboth{$ctg} $count{$ctg}\n";
  }
}
