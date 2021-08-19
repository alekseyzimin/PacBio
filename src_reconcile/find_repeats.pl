#!/usr/bin/env perl
#this code identifies repetitive contigs by coverage and multiple links
#coords on stdin
my $cov_thresh=2;
my $coords=$ARGV[0];
open(FILE,$coords);
#first compute coverage
while($line=<FILE>){
  chomp($line);
  push (@coords,$line);
  @F=split(/\s+/,$line);
  $cov{$F[-2]}+=($F[1]-$F[0]+1);
  $len{$F[-2]}=$F[11];
}

my $links=$ARGV[1];
open(FILE,$links);
while($line=<FILE>){
  my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$line);
  if($dir1 eq "F"){
    $edge_fwd{$ctg1}.="$ctg2 $dir2 $gap ";
    if($dir2 eq "F"){
      $edge_rev{$ctg2}.="$ctg1 F $gap ";
    }else{
      $edge_fwd{$ctg2}.="$ctg1 R $gap ";
    }
  }else{
    my $tdir=($dir2 eq "F") ? "R" : "F";
    $edge_rev{$ctg1}.="$ctg2 $tdir $gap ";
    if($dir2 eq "F"){
      $edge_rev{$ctg2}.="$ctg1 R $gap ";
    }else{
      $edge_fwd{$ctg2}.="$ctg1 F $gap ";
    }
  }
}

#identify repeats by high coverage and >2 in and >2 out
foreach $c(keys %edge_fwd){
  if(defined($edge_rev{$c})){
    my @f1=split(/\s+/,$edge_fwd{$c});
    my @f2=split(/\s+/,$edge_rev{$c});
    print "$c rev $edge_rev{$c}\n$c fwd $edge_fwd{$c}\n" if($cov{$c}/($len{$c}+1)>=$cov_thresh && $#f1>2 && $#f2>2);#+1 avoids /0 error
  }
}
