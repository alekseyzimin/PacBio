#!/usr/bin/env perl
#
#this code finds a set of non-overlapping matches for a mega-read
#
#
#assume input sorted by mega-read size

@interval_bgn=();
@interval_end=();
@interval_g_bgn=();
@interval_g_end=();

@mr=();
@nsr=();
@mrseq=();
@mrname=();
#max overlap is in percentage of size
$max_overlap_pct=$ARGV[0];
$kmer=$ARGV[1];
$fudge_factor=1.2;

while($line=<STDIN>){
    chomp($line);
($pbgn,$pend,$mbgn,$mend,$qlt,$scr,$mrlen,$pb,$mrseq,$mrname)=split(/\s+/,$line);
push(@intervals_sorted,"$pbgn $pend $mbgn $mend $mrlen $pb $mrseq $mrname");
}

#here we merge the matches
@merged_intervals=();
@curr_intervals=();
@f=split(/\s+/,$intervals_sorted[0]);
$last_mr=$f[6];
push(@curr_intervals,$intervals_sorted[0]);
for($i=1;$i<=$#intervals_sorted;$i++){
    @f=split(/\s+/,$intervals_sorted[$i]);
if(not($f[6] eq $last_mr)){
    if($#curr_intervals==0){
        push(@merged_intervals,$curr_intervals[0]);
    }else{
        $merge_index=0;
        $curr_intervals_output[$merge_index]=$curr_intervals[0];
        for($j=1;$j<=$#curr_intervals;$j++){
            @ff1=split(/\s+/,$curr_intervals_output[$merge_index]);
            @ff2=split(/\s+/,$curr_intervals[$j]);
            if(($ff2[0]-$ff1[1])-($ff2[2]-$ff1[3])>-.02*(($ff2[0]-$ff1[1])+($ff2[2]-$ff1[3])) && $ff2[2]-$ff1[3]>-5){ # if insertion in pb
                $curr_intervals_output[$merge_index]="$ff1[0] $ff2[1] $ff1[2] $ff2[3] $ff1[4] $ff1[5] $ff1[6] $ff1[7]";
            }else{
                $merge_index++;
                $curr_intervals_output[$merge_index]=$curr_intervals[$j];
            }
        }
    for($j=0;$j<=$#curr_intervals_output;$j++){
        push(@merged_intervals,$curr_intervals_output[$j]);
    }
    }
    @curr_intervals=();
}
$last_mr=$f[6];
push(@curr_intervals,$intervals_sorted[$i]);
}
    if($#curr_intervals==0){
        push(@merged_intervals,$curr_intervals[0]);
    }else{
        $merge_index=0;
        $curr_intervals_output[$merge_index]=$curr_intervals[0];
        for($j=1;$j<=$#curr_intervals;$j++){
            @ff1=split(/\s+/,$curr_intervals_output[$merge_index]);
            @ff2=split(/\s+/,$curr_intervals[$j]);
            if(($ff2[0]-$ff1[1])-($ff2[2]-$ff1[3])>-.02*(($ff2[0]-$ff1[1])+($ff2[2]-$ff1[3])) && $ff2[2]-$ff1[3]>-5){ # if insertion in pb
                $curr_intervals_output[$merge_index]="$ff1[0] $ff2[1] $ff1[2] $ff2[3] $ff1[4] $ff1[5] $ff1[6] $ff1[7]";
            }else{
                $merge_index++;
                $curr_intervals_output[$merge_index]=$curr_intervals[$j];
            }
        }
    for($j=0;$j<=$#curr_intervals_output;$j++){
        push(@merged_intervals,$curr_intervals_output[$j]);
    }
    }

#done merging, result is in @merged_intervals

foreach $interval(@merged_intervals){
($bgn,$end,$mbgn,$mend,$mrlen,$pb,$mrseq,$mrname)=split(/\s+/,$interval);
$max_overlap=$max_overlap_pct*($mend-$mbgn+1)/100;
$max_overlap=$kmer*$fudge_factor if($max_overlap<$kmer*$fudge_factor);
$overlap=0;
for($i=0;$i<=$#interval_g_bgn;$i++){

    last if($bgn >= $interval_g_bgn[$i] && $end <= $interval_g_end[$i]);#contained
    last if($bgn < $interval_g_bgn[$i] && $end > $interval_g_end[$i]);#containing -- happens on rare occasions

    $bgn_inside=0;
    $end_inside=0;
    $bgn_inside=1 if($bgn>=$interval_g_bgn[$i] && $bgn<=$interval_g_end[$i]);
    $end_inside=1 if($end>=$interval_g_bgn[$i] && $end<=$interval_g_end[$i]);
    next if($bgn_inside==0 && $end_inside==0);#non-overlapping this interlal
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
if($i>$#interval_g_bgn){#not contained anywhere
    if($overlap==0){
	push(@interval_g_bgn,$bgn);
	push(@interval_g_end,$end);
    }
    push(@intervals_out,"$bgn $end $mbgn $mend $mrlen $pb $mrseq $mrname");
}
}

@intervals_out_sorted= sort by_first_field @intervals_out;

foreach $interval(@intervals_out_sorted){
print $interval,"\n";
}

sub by_first_field 
{
    @f1=split(/\s+/,$a);
    @f2=split(/\s+/,$b);
    return($f1[0] <=> $f2[0]);
}
