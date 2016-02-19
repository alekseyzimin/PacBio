#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
$prev_match = "";
$match_ref_beg = 0;
$match_ref_end = 0;
$match_qry_beg = 0;
$match_qry_end = 0;
if(defined($ARGV[0])){
$max_gap_diff=$ARGV[0];
}else{
$max_gap_diff=500;
}

while ($line = <STDIN>) {
    chomp($line); 
    $line=~s/^\s+//;
    @currentFlds = split(" ",$line); 
    $currentMidpoint = ($currentFlds[3]+$currentFlds[4]) / 2;
    $local_direction = &reportDirectionBasedOnOrderedCoords (@currentFlds[3..4]);
    $keepMatchLine = 0;
    $current_match = "$currentFlds[17] $currentFlds[18]";
    if($current_match eq $prev_match) {
	next if ($local_direction != $match_direction);
	if (($prevFlds[3] < $prevFlds[4]) && ($prevMidpoint < $currentMidpoint)) {
            push(@matches,join(" ",@currentFlds));
	    my $gap_diff=abs(($currentFlds[0]-$prevFlds[1])-($currentFlds[3]-$prevFlds[4])); #this is the gap difference
	    my $covered = ($prevFlds[1]-$prevFlds[0] + $currentFlds[1] - $currentFlds[0]);
	    if($gap_diff >$max_gap_diff || 2*$covered < $currentFlds[0] - $prevFlds[1]){
                $badJoin=1;
            }else{
                $keepMatchLine = 1;
            }
	}elsif (($prevFlds[3] >= $prevFlds[4]) && ($prevMidpoint >= $currentMidpoint)) {
	    push(@matches,join(" ",@currentFlds));
            my $gap_diff = abs(($currentFlds[0]-$prevFlds[1])-($prevFlds[4]-$currentFlds[3])); #this is the gap difference
	    my $covered = ($prevFlds[1]-$prevFlds[0] + $currentFlds[1] - $currentFlds[0]);
            if($gap_diff >$max_gap_diff || 2*$covered < $currentFlds[0] - $prevFlds[1]){
	        $badJoin=1;
	    }else{
		$keepMatchLine = 1; 
	    }
	}
    }else{ # The ref or query has changed (or is the first line)
	if ($prev_match ne "") { # Output if not the first line
	if($badJoin){   
	   foreach $v(@matches){
	   print "$v\n";
	   }
	}else{
	    &outputMatchGroup;
	}
	}
	# Update the fields needed for only for a new match group (the others are updated below)
	$keepMatchLine = 1;
	$match_ref_beg = $currentFlds[0]; 
	$match_qry_beg = $currentFlds[3];
	$prev_match = $current_match;
	$match_direction = $local_direction;
	$match_bases = 0;
	$matching_bases=0;
	$badJoin=0;
	@matches=();
	push(@matches,join(" ",@currentFlds));
    }
    if ($keepMatchLine) {
	$match_ref_end = $currentFlds[1];
	$match_qry_end = $currentFlds[4];
	$matching_bases += $currentFlds[7]*$currentFlds[9]/100;
	$match_bases += $currentFlds[7];
	@prevFlds = @currentFlds;
	$prevMidpoint = $currentMidpoint;
    }
}

# Output for the last match (group)
if ($prev_match ne "") {
        if($badJoin){   
           foreach $v(@matches){
           print "$v\n";
           } 
        }else{
            &outputMatchGroup;
        }
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
    print "$match_ref_beg $match_ref_end | $match_qry_beg $match_qry_end | $ref_match_len $qry_match_len | $pctIdentityStr | @prevFlds[11..12] | $pctRefMatchLenStr $pctQueryMatchLenStr | @prevFlds[17..18]\n";
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


    
	
