#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $pbseqfile=$ARGV[0];
my $max_gap=$ARGV[1];
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
my $outread="";
my $last_coord =-1000000000;
#now we process the pb+mega-reads file
while($line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
	if(not($outread eq "")){
	    $indx=0;
	    @f=split(/(N)\1+/,$outread);
	    for($i=0;$i<=$#f;$i+=2){
	    	print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>100);
		$indx++;
	    }
	}
	$outread="";
	$last_coord =-1000000000;
	$rn=substr($line,1);

    }else{
	($bgn,$end,$pb,$seq,$name)=split(/\s+/,$line);
	if($outread eq ""){
		$outread=$seq;
	}else{
	    if($bgn>$last_coord){
		my $min_len=length($outread)<length($seq)?length($outread):length($seq);
		my $max_gap_local=int($min_len*.15);
		$max_gap_local=$max_gap if($max_gap_local>$max_gap);
		$max_gap_local=5 if($max_gap_local<5);
		if($bgn-$last_coord>$max_gap_local){#then put N's and later split
		$outread.=("N" x ($bgn-$last_coord));
		}else{
		$outread.=lc(substr($pbseq{$rn},$last_coord+1,$bgn-$last_coord));
		}
		$outread.=$seq;
	    }else{#overlapping
 	    $outread.=substr($seq,$last_coord-$bgn+1);
	    }
	}
    $last_coord=$end;
    }
}
#output the last one
            $indx=0;
            @f=split(/(N)\1+/,$outread);
            for($i=0;$i<=$#f;$i+=2){
                print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n";
                $indx++;
            }

