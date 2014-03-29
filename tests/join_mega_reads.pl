#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $pbseqfile=$ARGV[0];
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
	    $outread=substr($outread,0,length($pbseq{$rn})) if(length($outread)>length($pbseq{$rn})); #trim if longer than pb read
	    print ">$rn\n$outread\n";
	}
	$outread="";
	$last_coord =-1000000000;
	$rn=substr($line,1);

    }else{
	@f=split(/\s+/,$line);
	if($outread eq ""){
	    if($f[1]<0){#trim the front if sticking out of pb read
		$outread.=substr($f[4],-1*$f[1]);
	    }else{
		$outread.=$f[4];
	    }
	}else{
	    if($f[1]>$last_coord){
		$outread.=substr($pbseq{$rn},$last_coord,$f[1]-$last_coord);
	    }
	    $outread.=$f[4];
	}
	$last_coord=$f[2];
    }
}

