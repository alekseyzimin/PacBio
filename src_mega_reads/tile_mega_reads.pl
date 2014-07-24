#!/usr/bin/env perl
#
#this code finds a set of non-overlapping matches for a mega-read
#
#
#assume input sorted by mega-read size

use strict;
use warnings;

my @interval_bgn=();
my @interval_end=();
my @interval_g_bgn=();
my @interval_g_end=();

my @mr=();
my @nsr=();
my @mrseq=();
my @mrname=();
#max overlap is in percentage of size
my $max_overlap_pct=$ARGV[0];
my $kmer=$ARGV[1];
my $fudge_factor=1.2;

while(my $line=<STDIN>){
  chomp($line);
  my ($index, $pbgn, $pend, $lpath, $ldensity, $unitigs, $mrseq, $mrlen)=split(/\s+/,$line);
  my $max_overlap=$max_overlap_pct*length($mrseq)/100;
  $max_overlap=$kmer*$fudge_factor if($max_overlap<$kmer*$fudge_factor);

  my $bgn=$pbgn;
  my $end=$pend;
  my $overlap=0;
  my $i;
  for($i = 0; $i < @interval_g_bgn; $i++){
    last if($bgn >= $interval_g_bgn[$i] && $end <= $interval_g_end[$i]);#contained
    last if($bgn < $interval_g_bgn[$i] && $end > $interval_g_end[$i]);#containing -- happens on rare occasions

    my $bgn_inside=0;
    my $end_inside=0;
    $bgn_inside=1 if($bgn>=$interval_g_bgn[$i] && $bgn<=$interval_g_end[$i]);
    $end_inside=1 if($end>=$interval_g_bgn[$i] && $end<=$interval_g_end[$i]);
    next if($bgn_inside==0 && $end_inside==0);#non-overlapping this interval
#we get here if we have an overlap
    if($bgn_inside==1){#check the overlaps
      last if($interval_g_end[$i]-$bgn>$max_overlap);#overlap bigger -- contained
      $interval_g_end[$i]=$end;
      $overlap=1;
    }else{
      last if($end-$interval_g_bgn[$i]>$max_overlap);
      $interval_g_bgn[$i]=$bgn;
      $overlap=1;
    }
  }
#  print("mo:$max_overlap ($pbgn, $pend) ovl:$overlap ", $i < @interval_g_bgn ? "($interval_g_bgn[$i],$interval_g_end[$i])" : "-", "\n");

  if($i>$#interval_g_bgn){#not contained anywhere
    if($overlap==0){
      push(@interval_g_bgn,$bgn);
      push(@interval_g_end,$end);
    }
    push(@interval_bgn,$bgn);
    push(@interval_end,$end);
    push(@nsr,$lpath);
    # push(@mrseq,$mrseq);
    # push(@mrnames,$mrname);
  }
}

if(@interval_bgn) {
  for(my $i=0; $i < @interval_bgn; $i++) {
    print("$interval_bgn[$i] $interval_end[$i] $nsr[$i]\n")
  }
}
