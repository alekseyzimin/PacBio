#!/usr/bin/env perl
#
$nmatches=2;
$score_window=.95;
$min_overlap=1;
$start=0;
$end=0;
@current_matches=();
@current_ends=();
@bestmtch=(); 
while($line=<STDIN>){
    chomp($line);
    @f=split(/\s+/,$line);
    if($f[0]>$start){
        if(scalar(@bestmtch)>0){
	    @ff=split(/\s+/,$bestmtch[$#bestmtch]);
	    $max_score=$ff[8]*$score_window;
	    $bgn=$#bestmtch-$nmatches>=0?$#bestmtch-$nmatches:0;
	    for($i=$bgn;$i<=$#bestmtch;$i++){
		@ff=split(/\s+/,$bestmtch[$i]);
		if($ff[8]>=$max_score){
		    $end=$ff[1];
		    print $bestmtch[$i],"\n";
		    push(@current_matches,$ff[12]);
		    push(@current_ends,$ff[1]);
		}
	    }
	    @bestmtch=();
        }
        $start=$f[0];
    }
    push(@bestmtch,$line) if(check_match($f[0],$f[1],$end,$f[12]));   
}


sub check_match{
    my ($first_coord,$last_coord,$last_end,$name)=@_;
    my $flag=0;
    if($last_coord>=$last_end){
            $flag=1;
	    for(my $i=0;$i<=$#current_matches;$i++){
		if($current_ends[$i]>$first_coord+$min_overlap){
		    $flag=0;
		    if(overlap($current_matches[$i],$name)){
			$flag=1;
			last;
		    }
		}
	    }
        }
    return($flag);
}


sub overlap{
    my @f1=split(/_/,$_[1]);
    my $flag=0;
    for(my $j=$#f1;$j>=0;$j--){
	$tmpstr=join("_",@f1[0..$j]);
	if(index($_[0],$tmpstr)>0 && index($_[0],$tmpstr)==length($_[0])-length($tmpstr)){
	    $flag=1;
	    last;
	}
    }
    return($flag);
}
