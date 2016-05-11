#!/usr/bin/env perl
my $match_ref_beg = 0;
my $match_ref_end = 0;
my $match_qry_beg = 0;
my $match_qry_end = 0;
my @output_matches = ();
my $max_gap_diff =500;
if(defined($ARGV[0])){
$max_gap_diff=$ARGV[0];
}

my %ctg_lines=();
my $scf="";
while($line=<STDIN>){
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  if(not($f[-2] eq $scf)){
    if(not($scf eq "")){
      my $ctg;
      foreach $ctg( keys %ctg_lines ){
      merge_matches(split("\n",$ctg_lines{$ctg}));
      }
      print sort by_first_field @output_matches;
      @output_matches = ();
    }
    %ctg_lines=();
    $scf=$f[-2];
    }
  $ctg_lines{$f[-1]}.=$line."\n";
}
my $ctg;
foreach $ctg( keys %ctg_lines ){
  merge_matches(split("\n",$ctg_lines{$ctg}));
  }
@sorted= sort by_first_field @output_matches;
print @sorted;



sub merge_matches{

my $prev_match = "";
my @currentFlds =();
my $currentMidpoint;
my $prevMidpoint;
my $local_direction;
my $match_direction = 0;
my $keepMatchLine = 0;
my $groupNotEmpty = 0;
#assumes matches for one contig at a time
foreach my $line (@_) {
    @currentFlds = split(/\s+/,$line); 
    $currentMidpoint = ($currentFlds[3]+$currentFlds[4]) / 2;
    if($keepMatchLine == 1){#we check the next match
      $local_direction = &reportDirectionBasedOnOrderedCoords (@currentFlds[3..4]);
      $keepMatchLine = 0;
      if ($local_direction == $match_direction){
        if (($prevFlds[3] < $prevFlds[4]) && ($prevMidpoint < $currentMidpoint)) { 
          $keepMatchLine = 1 if(abs(($currentFlds[0]-$prevFlds[1])-($currentFlds[3]-$prevFlds[4])) <= $max_gap_diff);
        }elsif (($prevFlds[3] >= $prevFlds[4]) && ($prevMidpoint >= $currentMidpoint)) { 
          $keepMatchLine = 1 if(abs(($currentFlds[0]-$prevFlds[1])-($prevFlds[4]-$currentFlds[3])) <= $max_gap_diff);
        }
      }
    }
    
    if($keepMatchLine == 0){
      &outputMatchGroup if($groupNotEmpty);
      $match_direction = &reportDirectionBasedOnOrderedCoords (@currentFlds[3..4]);
      $match_ref_beg = $currentFlds[0];
      $match_qry_beg = $currentFlds[3];
      $match_bases = $currentFlds[7];
      $matching_bases = $currentFlds[7]*$currentFlds[9]/100;
      $keepMatchLine = 1;
      $groupNotEmpty = 1;
    }else{
      $matching_bases += $currentFlds[7]*$currentFlds[9]/100;
      $match_bases += $currentFlds[7];
    }
    $match_ref_end = $currentFlds[1];
    $match_qry_end = $currentFlds[4];
    @prevFlds = @currentFlds;
    $prevMidpoint = $currentMidpoint;
}
&outputMatchGroup if($groupNotEmpty);
}

sub reportDirectionBasedOnOrderedCoords
{
    my ($val1, $val2) = @_;
    my ($result);

    if ($val1 < $val2) {
	$result = 1; }
    else {
	$result = -1; }

    return ($result);
}

sub outputMatchGroup
{
    $qry_match_len = abs($match_qry_end-$match_qry_beg) + 1;
    $ref_match_len = $match_ref_end-$match_ref_beg + 1;
    $pctIdentity = $matching_bases*100 / $match_bases;
    $pctRefMatchLen = 100 * ($ref_match_len/$prevFlds[11]);
    $pctQueryMatchLen = 100 * ($qry_match_len/$prevFlds[12]);
    $pctIdentityStr = &makeHundredths ($pctIdentity);
    $pctRefMatchLenStr = &makeHundredths ($pctRefMatchLen);
$pctQueryMatchLenStr = &makeHundredths ($pctQueryMatchLen);
    push(@output_matches,"$match_ref_beg $match_ref_end | $match_qry_beg $match_qry_end | $ref_match_len $qry_match_len | $pctIdentityStr | @prevFlds[11..12] | $pctRefMatchLenStr $pctQueryMatchLenStr | @prevFlds[17..18]\n");
}

sub makeHundredths
{
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
    
	
