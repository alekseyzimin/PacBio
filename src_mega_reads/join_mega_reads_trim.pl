#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $pbseqfile=$ARGV[0];
my $max_gap=$ARGV[1];
my $allowed_gaps=$ARGV[2];
my $kmer=$ARGV[3];
my $good_pb=$ARGV[4];
my $fudge_factor=1.2;
$kmer*=$fudge_factor;

my $rn="";
my %pbseq;
open(FILE,$pbseqfile);
while($line=<FILE>){
    chomp($line);
    if($line =~ /^>/){
	@f=split(/\s+/,$line);
	$rn=substr($f[0],1);
    }else{
	$pbseq{$rn}.=$line;
    }
}

open(FILE,$allowed_gaps);
while($line=<FILE>){
    chomp($line);
    @f=split(/\s+/,$line);
    $allowed{"$f[0] $f[1]"}=1;
}

open(FILE,$good_pb);
while($line=<FILE>){
    chomp($line);
    $good_pb{$line}=1;
}

my @lines=();
my $outread="";
#now we process the pb+mega-reads file
while($line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
	if(@lines){
	    @lines_sorted = sort by_first_number @lines;
	    $outread=&process_sorted_lines(@lines_sorted);
	    if(not($outread eq "")){
		$indx=0;
		@f=split(/N/,$outread);
		for($i=0;$i<=$#f;$i++){
		    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=400);
		    $indx++;
		}
	    }
	}
	@lines=();
    }else{
	push(@lines, $line);
    }
}
#do not forgrt the last one
        if(@lines){
            @lines_sorted = sort by_first_number @lines;
            $outread=&process_sorted_lines(@lines_sorted);
            if(not($outread eq "")){
                $indx=0;
                @f=split(/N/,$outread);
                for($i=0;$i<=$#f;$i++){
                    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=400);
                    $indx++;
                }
            }
        }



sub by_first_number{
($fn1,$rest)=split(/\s+/,$a);
($fn2,$rest)=split(/\s+/,$b);
return($a <=> $b);
}

sub process_sorted_lines{
    my $outread="";
    $last_coord =-1000000000;
    foreach $l(@_){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=split(/\s+/,$l);
        $seq=substr($mseq,$mbgn-1,$mend-$mbgn+1);
        die("inconsistent sequence length") if(not(length($mseq)==$mlen));
        die("pacbio read $pb does not exist in the sequence file!!!") if(not(defined($pbseq{$rn})));

        if($outread eq ""){
	    $outread=$seq;
        }else{
            my @k1s=split(/_/,$last_mr);
            my @k2s=split(/_/,$name);
            $join_allowed=0;
            $k1s[$#k1s] =substr($k1s[$#k1s],0,-1);
            $k2s[0] = substr($k2s[0],0,-1);
            $str="$k1s[$#k1s] $k2s[0]";
            $str="$k2s[0] $k1s[$#k1s]" if($k1s[$#k1s]>$k2s[0]);
            $join_allowed=1 if($allowed{$str});
            $join_allowed=1 if(defined($good_pb{$pb}));

            if($bgn>$last_coord){#if gap -- check if the closure is allowed
                my $min_len=0;
                my $max_gap_local;

                if(defined($good_pb{$pb})){
		    $max_gap_local=length($outread)>length($seq)?length($outread):length($seq);
                }else{
		    $min_len=length($outread)<length($seq)?length($outread):length($seq);
		    $max_gap_local=int($min_len*0.3);
		    $max_gap_local=$max_gap if($max_gap_local>$max_gap);
		    $max_gap_local=25 if($max_gap_local<25);
                }
                if($bgn-$last_coord<$max_gap_local && $join_allowed){#then put N's and later split
		    $outread.=lc(substr($pbseq{$rn},$last_coord+1,$bgn-$last_coord)).$seq;
                }else{
		    $outread.="N".$seq;
                }
            }else{#overlapping
		$join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads 
		$join_allowed=1 if($last_implied_coord-($bgn-$mbgn+1)>1 && $last_implied_coord-($bgn-$mbgn+1)<=$kmer);
		if($join_allowed){
		    $outread.=substr($seq,$last_coord-$bgn+1);
		}else{
		    $outread.="N".$seq;
		}
            }
        }
	$last_coord=$end;
	$last_implied_coord=$end+$mlen-$mend;
	$last_mr=$name;
	$last_len=length($seq);
	$last_ext=substr($mseq,$mend);
	$last_mend=$mend;
    }
    return($outread);
}
