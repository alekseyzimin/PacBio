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
	    $outread=substr($outread,0,length($outread)-($last_coord-length($pbseq{$rn}))) if($last_coord>length($pbseq{$rn})); #trim if longer than pb read
	    $indx=0;
	    @f=split(/(N)\1+/,$outread);
	    for($i=0;$i<=$#f;$i+=2){
	    	print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n";
		$indx++;
	    }
	}
	$outread="";
	$last_coord =-1000000000;
	$rn=substr($line,1);

    }else{
	($bgn,$end,$pb,$seq)=split(/\s+/,$line);
	if($outread eq ""){
	    if($bgn<0){#trim the front if sticking out of pb read
		$outread=substr($,-1*$bgn);
	    }else{
		$outread=$seq;
	    }
	}else{
	    if($bgn>$last_coord){
		if($bgn-$last_coord>$max_gap){#then put N's and later split
		$outread.=("N" x ($f[0]-$last_coord));
		}else{
		$outread.=lc(substr($pbseq{$rn},$last_coord+1,$f[0]-$last_coord-1));
		}
		$outread.=$f[4];
	    }else{
 	    $outread.=substr($f[4],$last_coord-$f[1]+1);
	    #$outread.=$f[4];
	    }
	}
	$last_coord=$end;
    }
}
#output the last one
           $outread=substr($outread,0,length($pbseq{$rn})) if(length($outread)>length($pbseq{$rn})); #trim if longer than pb read
            $indx=0;
            @f=split(/(N)\1+/,$outread);
            for($i=0;$i<=$#f;$i+=2){
                print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n";
                $indx++;
            }

