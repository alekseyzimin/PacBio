#!/usr/bin/env perl
#
#this code finds a set of non-overlapping matches for a mega-read
#
#
#assume input sorted by mega-read size
#max overlap is in percentage of size
$max_overlap_pct=$ARGV[0];
$kmer=$ARGV[1];
$seqfile=$ARGV[2];
$mr_namefile=$ARGV[3];
$min_len=0;
$min_len=$ARGV[4] if(not($ARGV[4] eq ""));

open(FILE,$seqfile);
while($line=<FILE>){
    chomp($line);
    if($line =~ /^\>/){
	$rn=substr($line,1);
    }else{
	$seq{$rn}=$line;
	$line=reverse($line);
	$line=~tr/ACGTNacgtn/TGCANtgcan/;
	$seq{$rn+1}=$line;
    }
}

@mr_ext_name=();
open(FILE,$mr_namefile);
while($line=<FILE>){
    chomp($line);
    push(@mr_ext_name,$line);
}

$last_pb_read="";
@lines=();
while($l=<STDIN>){
    chomp($l);
    my @ff=split(/\s+/,$l);
    next if($ff[10] eq "");
    next if(($ff[7]-$ff[6])<$min_len);
    #next if($ff[3]==1);
    $original_mega_read=$ff[1];
    $mega_read=$original_mega_read;
    die("mega-read $mega_read from $ff[1] has no sequence!") if(not(defined($seq{$mega_read})));
    @fff=split(/\//,$ff[0]);
    $pb_read=join("/",@fff[0..($#fff-1)]);
    $pb_len{$pb_read}=$ff[11];
    $sequence=$seq{$mega_read};
    if($ff[3]==1){
        @kunis=split(/_/,$mega_read);
        $mega_read=join("_",reverse(@kunis));
        $mega_read=~tr/FR/RF/;
        $sequence=reverse($sequence);
        $sequence=~tr/ACGTNacgtn/TGCANtgcan/;
    }
    if(not($pb_read eq $last_pb_read)){
	if(not($last_pb_read eq "")){
	    print ">$last_pb_read\n";
	    create_tiling(sort {$$a[7] <=> $$b[7] || $$a[0] <=> $$b[0]} @lines) if(@lines);
	}
	@lines=();
	$last_pb_read=$pb_read;
    }
    $mtch_bases=($ff[7]-$ff[6])*$ff[5]/100;
    $weight=($ff[7]-$ff[6])/(101-$ff[5]);
    #print "DEBUG $original_pb $ff[0] $ff[9] $ff[10] $ff[6] $ff[7] $ff[1]\n";
    push(@lines,[$ff[9],$ff[10],$ff[6],$ff[7],$ff[8],$pb_read,$sequence,$mega_read,$mtch_bases,$weight]);
}
#last one
print ">$last_pb_read\n";
create_tiling(sort {$$a[7] <=> $$b[7] || $$a[0] <=> $$b[0]} @lines) if(@lines);

sub merge_intervals_detect{
    my @intervals_sorted=@_;
    my @merged_intervals=();
    my @curr_intervals=();
    my @curr_intervals_output=();
    my @f=@{$intervals_sorted[0]};
    my $last_mr=$f[7];
    push(@curr_intervals,$intervals_sorted[0]);
    for($i=1;$i<=$#intervals_sorted;$i++){
        @f=@{$intervals_sorted[$i]};
        if(not($f[7] eq $last_mr)){
            $merge_index=0;
            push(@curr_intervals_output,$curr_intervals[0]);
            for($j=1;$j<=$#curr_intervals;$j++){
                @ff1=@{$curr_intervals_output[$merge_index]};
                @ff2=@{$curr_intervals[$j]};
                my $covered=$ff1[1]-$ff1[0]+$ff2[1]-$ff2[0];
                my $gap_pb=$ff2[0]-$ff1[1];
                my $gap_mr=$ff2[2]-$ff1[3];
                if(abs($gap_pb-$gap_mr)<5000 && 2*$covered>=$gap_pb && $ff2[2]>$ff1[2]){ #if insertion in pb
                    my $qlt=$ff1[8]+$ff2[8];
                    $curr_intervals_output[$merge_index]=[$ff1[0],$ff2[1],$ff1[2],$ff2[3],$ff1[4..7],$qlt,$ff1[9]]; #do the merging
                }else{
                    $merge_index++;
                    $curr_intervals_output[$merge_index]=$curr_intervals[$j];
                }
            }
            foreach $interval(@curr_intervals_output){
                push(@merged_intervals,$interval);
            }
            $last_mr=$f[7];
            @curr_intervals=();
            @curr_intervals_output=();
        }
        push(@curr_intervals,$intervals_sorted[$i]);
    }
    $merge_index=0;
    push(@curr_intervals_output,$curr_intervals[0]);
    for($j=1;$j<=$#curr_intervals;$j++){
	@ff1=@{$curr_intervals_output[$merge_index]};
	@ff2=@{$curr_intervals[$j]};
	my $covered=$ff1[1]-$ff1[0]+$ff2[1]-$ff2[0];
	my $gap_pb=$ff2[0]-$ff1[1];
	my $gap_mr=$ff2[2]-$ff1[3];
	if(abs($gap_pb-$gap_mr)<5000 && 2*$covered>=$gap_pb && $ff2[2]>$ff1[2]){ #if insertion in pb
	    my $qlt=$ff1[8]+$ff2[8];
	    $curr_intervals_output[$merge_index]=[$ff1[0],$ff2[1],$ff1[2],$ff2[3],$ff1[4..7],$qlt,$ff1[9]]; #do the merging
	}else{
	    $merge_index++;
	    $curr_intervals_output[$merge_index]=$curr_intervals[$j];
	}
    }
    foreach $interval(@curr_intervals_output){
	push(@merged_intervals,$interval);
    }
    
    return(@merged_intervals);
}


sub merge_intervals{
    my @intervals_sorted=@_;
    my @merged_intervals=();
    my @curr_intervals=();
    my @curr_intervals_output=();
    my @f=@{$intervals_sorted[0]};
    my $last_mr=$f[7];
    push(@curr_intervals,$intervals_sorted[0]);
    for($i=1;$i<=$#intervals_sorted;$i++){
        @f=@{$intervals_sorted[$i]};
        if(not($f[7] eq $last_mr)){
	    $merge_index=0;
	    push(@curr_intervals_output,$curr_intervals[0]);
	    for($j=1;$j<=$#curr_intervals;$j++){
		@ff1=@{$curr_intervals_output[$merge_index]};
		@ff2=@{$curr_intervals[$j]};
		my $covered=$ff1[1]-$ff1[0]+$ff2[1]-$ff2[0];
		my $gap_pb=$ff2[0]-$ff1[1];
		my $gap_mr=$ff2[2]-$ff1[3];
		if(abs($gap_pb-$gap_mr)<5000 && ($gap_mr <= $gap_pb-5 || $gap_mr <= $gap_pb*1.1) && $gap_pb>-5 && 2*$covered>=$gap_pb && $ff2[2]>$ff1[2]){ #if insertion in pb
		    my $qlt=$ff1[8]+$ff2[8];
		    $curr_intervals_output[$merge_index]=[$ff1[0],$ff2[1],$ff1[2],$ff2[3],$ff1[4..7],$qlt,$ff1[9]]; #do the merging
		}else{
		    $merge_index++;
		    $curr_intervals_output[$merge_index]=$curr_intervals[$j];
		}
	    }
	    foreach $interval(@curr_intervals_output){
		push(@merged_intervals,$interval);
	    }
	    $last_mr=$f[7];
	    @curr_intervals=();
	    @curr_intervals_output=();
        }
        push(@curr_intervals,$intervals_sorted[$i]);
    }
    $merge_index=0;
    push(@curr_intervals_output,$curr_intervals[0]);
    for($j=1;$j<=$#curr_intervals;$j++){
	@ff1=@{$curr_intervals_output[$merge_index]};
	@ff2=@{$curr_intervals[$j]};
	my $covered=$ff1[1]-$ff1[0]+$ff2[1]-$ff2[0];
	my $gap_pb=$ff2[0]-$ff1[1];
	my $gap_mr=$ff2[2]-$ff1[3];
	if(abs($gap_pb-$gap_mr)<5000 && ($gap_mr <= $gap_pb-5 || $gap_mr <= $gap_pb*1.1) && $gap_pb>-5 && 2*$covered>=$gap_pb && $ff2[2]>$ff1[2]){ #if insertion in pb
	    my $qlt=$ff1[8]+$ff2[8];
	    $curr_intervals_output[$merge_index]=[$ff1[0],$ff2[1],$ff1[2],$ff2[3],$ff1[4..7],$qlt,$ff1[9]]; #do the merging
	}else{
	    $merge_index++;
	    $curr_intervals_output[$merge_index]=$curr_intervals[$j];
	}
    }
    foreach $interval(@curr_intervals_output){
	push(@merged_intervals,$interval);
    }

    return(@merged_intervals);
}

sub create_tiling{
    @interval_bgn=();
    @interval_end=();
    @interval_g_bgn=();
    @interval_g_end=();
    @intervals_out=();
    $fudge_factor=1.2;
    my @lines=@_;
#we output merged intervals for detection of bad pb
    foreach $interval(merge_intervals_detect(@lines)){
	($bgn,$end,$mbgn,$mend,$mrlen,$pb,$mrseq,$mrname,$qlt,$scr)=@{$interval};
        #print STDERR "LINE $bgn $end $mbgn $mend $mrlen\n";
	print STDERR "$pb $mrname 0 0 0 ",$qlt/($mend-$mbgn)*100," $mbgn $mend $mrlen $bgn $end $pb_len{$pb} 0\n";
    }

    @merged_intervals_sorted=sort {$$b[8] <=> $$a[8]} merge_intervals(@lines);
    foreach $interval(@merged_intervals_sorted){
($bgn,$end,$mbgn,$mend,$mrlen,$pb,$mrseq,$mrname,$qlt,$scr)=@{$interval};
$max_overlap=$max_overlap_pct*($mend-$mbgn+1)/100;
$max_overlap=$kmer*$fudge_factor if($max_overlap<$kmer*$fudge_factor);
$overlap=0;
#print "DEBUG MERGING $bgn $end $mbgn $mend $mrlen $pb $qlt,$scr\n";
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
#print "$i $#interval_g_bgn $overlap\n";
if($i>$#interval_g_bgn){#not contained anywhere
    if($overlap==0){
	push(@interval_g_bgn,$bgn);
	push(@interval_g_end,$end);
    }
    push(@intervals_out,[$bgn,$end,$mbgn,$mend,$mrlen,$pb,$mrseq,$mr_ext_name[$mrname]]);
}
    }

@intervals_out_sorted= sort {$$a[0] <=> $$b[0]} @intervals_out;

foreach $interval(@intervals_out_sorted){
    print join(" ",@{$interval}),"\n";
}
}

