#!/usr/bin/env perl
my $max_gap_diff =500;
my $max_gap_allowed=10000000;
if(defined($ARGV[0])){
  $max_gap_diff=$ARGV[0];
}
if(defined($ARGV[1])){
  $max_gap_allowed=$ARGV[1];
}


my %ctg_lines=();
my $scf="";
while($line=<STDIN>){
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  if(not($f[-2] eq $scf)){
    if(not($scf eq "")){
      foreach my $ctg( keys %ctg_lines ){
        merge_matches(split("\n",$ctg_lines{$ctg}));
      }
    }
    %ctg_lines=();
    $scf=$f[-2];
  }
  $ctg_lines{$f[-1]}.=$line;
}
my $ctg;
foreach $ctg( keys %ctg_lines ){
  merge_matches(split("\n",$ctg_lines{$ctg}));
}

sub merge_matches{
#this is an almost proper LCS code
#from the input we have to decide whether LCS is forward or reverse and output only one
my @matches_fwd=();
my @matches_rev=();
my $rname;
my $qname;
my $rlen;
my $qlen;
if(scalar(@_)==1){
  print $_[0],"\n"
}else{
  foreach my $line (@_){
    @currentFlds = split(/\s+/,$line);
    if($rname eq ""){
      $rname=$currentFlds[-2];
      $qname=$currentFlds[-1];
      $rlen=$currentFlds[11];
      $qlen=$currentFlds[12];
    }
    push(@matches_fwd,$line) if($currentFlds[3]<$currentFlds[4]);
    push(@matches_rev,$line) if($currentFlds[3]>$currentFlds[4]);
  }
  #now we find LCS for the forward matches and LCS for the reverse matches and figure out which one is longer
  #first forward
  my @fwd_rstarts=();
  my @fwd_rends=();
  my @fwd_qstarts=();
  my @fwd_qends=();
  my @fwd_lens=();
  my @fwd_quals=();
  my $fwd_len=0;
  if(scalar(@matches_fwd)>0){
  for($i=0;$i<$#matches_fwd;$i++){
    @line_i=split(/\s+/,$matches_fwd[$i]);
    $fwd_len+=$line_i[7]*$line_i[9]/100;
    push(@fwd_rstarts,$line_i[0]);
    push(@fwd_rends,$line_i[1]);
    push(@fwd_qstarts,$line_i[3]);
    push(@fwd_qends,$line_i[4]);
    push(@fwd_lens,$line_i[7]);
    push(@fwd_quals,$line_i[7]*$line_i[9]);
    for($j=$i+1;$j<=$#matches_fwd;$j++){
      @line_j=split(/\s+/,$matches_fwd[$j]);
      my $diff=abs($line_j[0]-$line_i[1]-$line_j[3]+$line_i[4]);
      if($diff>$max_gap_diff && $line_j[0]-$line_i[1]<$max_gap_allowed && $line_j[3]-$line_i[4]<$max_gap_allowed){
        $i=$j-1;
        $j=$#matches_fwd+2;
      }else{
        $fwd_rends[-1]=$line_j[1];
        $fwd_qends[-1]=$line_j[4];
        $fwd_lens[-1]+=$line_j[7];
        $fwd_quals[-1]+=$line_j[7]*$line_j[9];
        $fwd_len+=$line_j[7]*$line_j[9]/100;
      }
    $i=$#matches_fwd if($j==$#matches_fwd);
    }
  }
  }
  my @rev_rstarts=();
  my @rev_rends=();
  my @rev_qstarts=();
  my @rev_qends=();
  my @rev_lens=();
  my @rev_quals=();
  my $rev_len=0;
  if(scalar(@matches_rev)>0){
  for(my $i=0;$i<$#matches_rev;$i++){
    @line_i=split(/\s+/,$matches_rev[$i]);
    $rev_len=$line_i[7]*$line_i[9]/100;
    push(@rev_rstarts,$line_i[0]);
    push(@rev_rends,$line_i[1]);
    push(@rev_qstarts,$line_i[4]);
    push(@rev_qends,$line_i[3]);
    push(@rev_lens,$line_i[7]);
    push(@rev_quals,$line_i[7]*$line_i[9]);
    for(my $j=$i+1;$j<=$#matches_rev;$j++){
      @line_j=split(/\s+/,$matches_rev[$j]);
      my $diff=abs($line_j[0]-$line_i[1]-$line_i[4]+$line_j[3]);
      if($diff>$max_gap_diff && $line_j[0]-$line_i[1]<$max_gap_allowed && $line_i[4]-$line_j[3]<$max_gap_allowed){
        $i=$j-1;
        $j=$#matches_rev+1;
      }else{
        $rev_rends[-1]=$line_j[1];
        $rev_qends[-1]=$line_j[3];
        $rev_lens[-1]+=$line_j[7];
        $rev_quals[-1]+=$line_j[7]*$line_j[9];
        $rev_len+=$line_j[7]*$line_j[9]/100;
      }
    $i=$#matches_rev if($j==$#matches_rev);
    }
  }
  }
  if($fwd_len>$rev_len){
#match is forward
    for(my $i=0;$i<=$#fwd_rstarts;$i++){
      print "$fwd_rstarts[$i] $fwd_rends[$i] | $fwd_qstarts[$i] $fwd_qends[$i] | ",$fwd_rends[$i]-$fwd_rstarts[$i]," ",$fwd_qends[$i]-$fwd_qstarts[$i]," | ",makeHundredths($fwd_quals[$i]/$fwd_lens[$i])," | $rlen $qlen | ",makeHundredths($fwd_lens[$i]/$rlen*100)," ",makeHundredths($fwd_lens[$i]/$qlen*100)," | $rname $qname\n";
    }
  }else{
    for(my $i=0;$i<=$#rev_rstarts;$i++){
      print "$rev_rstarts[$i] $rev_rends[$i] | $rev_qstarts[$i] $rev_qends[$i] | ",$rev_rends[$i]-$rev_rstarts[$i]," ",-$rev_qends[$i]+$rev_qstarts[$i]," | ",makeHundredths($rev_quals[$i]/$rev_lens[$i])," | $rlen $qlen | ",makeHundredths($rev_lens[$i]/$rlen*100)," ",makeHundredths($rev_lens[$i]/$qlen*100)," | $rname $qname\n";
    }
  }
}
}

sub makeHundredths{
    my ($value) = @_;
    $value *= 100;
    $value = int ($value+.50001);
    while (length ($value) < 3) {
	$value = "0$value"; }
    substr ($value, -2, 0) = ".";
    return ($value);
}

sub by_first_field{
my @f1=split(/\s+/,$a);
my @f2=split(/\s+/,$b);
return($f1[0] <=> $f2[0]);
}
    
	
