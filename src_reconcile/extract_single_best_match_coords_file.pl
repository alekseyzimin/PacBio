#!/usr/bin/env perl
#
#this code extracts a single best query match from coords file -- the match is decided by the number of matching bases * match identity
#
my %best_matches=(), %best_scores=();
my $linenumber=0;
while($line=<STDIN>){
  chomp($line);
  push(@lines,$line);
  my @f=split(/\s+/,$line);
  my $score=$f[7]*$f[9];
  if(not(defined($best_matches{$f[-1]})) || (defined($best_matches{$f[-1]}) && $best_scores{$f[-1]}<$score)){
    $best_matches{$f[-1]}=$linenumber;
    $best_scores{$f[-1]}=$score;
  }
  $linenumber++;
}

$linenumber=0;
foreach my $l(@lines){
  my @f=split(/\s+/,$l);  
  print $l,"\n" if($best_matches{$f[-1]}==$linenumber);
  $linenumber++;
}
