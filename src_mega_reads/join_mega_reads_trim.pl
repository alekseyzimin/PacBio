#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $pbseqfile=$ARGV[0];
my $max_gap=$ARGV[1];
my $allowed_gaps=$ARGV[2];
my $kmer=$ARGV[3];

my $rn="";
my %pbseq;
open(FILE,$pbseqfile);
while($line=<FILE>){
    chomp($line);
    if($line =~ /^>/){
	$rn=substr($line,1);
    }else{
	$pbseq{$rn}.=$line;
    }
}

open(FILE,$allowed_gaps);
while($line=<FILE>){
    chomp($line);
    $allowed{$line}=1;
}

my $outread="";
my $last_coord =-1000000000;
#now we process the pb+mega-reads file
while($line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
	if(not($outread eq "")){
	    $indx=0;
	    @f=split(/N/,$outread);
	    for($i=0;$i<=$#f;$i++){
	    	print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>400);
		$indx++;
	    }
	}
	$outread="";
	$last_coord =-1000000000;
	$last_mr_len=0;
	$rn=substr($line,1);

    }else{
	($bgn,$end,$pb,$seq,$name)=split(/\s+/,$line);
	if($outread eq ""){
		$outread=$seq;
	}else{
	    #first we figure out if we can allow this join.  
            my @k1s=split(/_/,$last_mr);
            my @k2s=split(/_/,$name);
	    $join_allowed=0;
	    $k1s[$#k1s] =substr($k1s[$#k1s],0,-1);
            $k2s[0] = substr($k2s[0],0,-1);
            $str="$k1s[$#k1s] $k2s[0]";
            $str="$k2s[0] $k1s[$#k1s]" if($k1s[$#k1s]>$k2s[0]);
            $join_allowed=1 if($allowed{$str});

	    if($bgn>$last_coord){#if gap -- check if the closure is allowed
		my $min_len=length($outread)<length($seq)?length($outread):length($seq);
		$max_gap_local=int($min_len*.33);
		$max_gap_local=$max_gap if($max_gap_local>$max_gap);
		$max_gap_local=50 if($max_gap_local<50);
		#print "join status ",$bgn-$last_coord," $max_gap_local $join_allowed\n";
		if($bgn-$last_coord<$max_gap_local && $join_allowed){#then put N's and later split
		$outread.=lc(substr($pbseq{$rn},$last_coord+1,$bgn-$last_coord)).$seq;
		}else{
		$outread.="N".$seq;
		}
	    }else{#overlapping
 	    $join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads 
	    $join_allowed=1 if($last_coord-$bgn>=10 && $last_coord-$bgn<=2*$kmer);
            # we join if same mega-read, just fractured, or the overlap is less than kmer length, or join is allowed
 	    if($join_allowed){
	    $outread.=substr($seq,$last_coord-$bgn+1);
            }else{
            $outread.="N".$seq;
            }
	    }
	}
    $last_coord=$end;
    $last_mr_len=length($seq);
    $last_mr=$name;
    }
}
#output the last one
            $indx=0;
            @f=split(/(N)\1+/,$outread);
            for($i=0;$i<=$#f;$i+=2){
                print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n";
                $indx++;
            }

