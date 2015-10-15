#!/usr/bin/env perl
#
#this code refines the alignments of mega reads to pacbio reads using nucmer
use mummer;

#first we read in PB sequences
my $pbseqfile=$ARGV[0];
my $PREFIX=$ARGV[1];
my %needed_pb;
my $rn="";
my %pbseq;
my $readnumber=0;
my %readnames;

while($line=<STDIN>){
    push(@file,$line);
    chomp($line);
    if($line =~ /^>/){
	($rn,$junk)=split(/\s+/,substr($line,1));
        $needed_pb{$rn}=1;
	}
}

open(FILE,$pbseqfile);
while($line=<FILE>){
    chomp($line);
    if($line =~ /^>/){
	@f=split(/\s+/,$line);
	$rn=substr($f[0],1);
    }else{
	$pbseq{$rn}.=$line if($needed_pb{$rn});
    }
}

print "pacbio mega-reads\nNUCMER\n";
my @lines=();
my $outread="";
open(OUTFILE1,">t.$PREFIX.maximal_mr.fa");
open(OUTFILE2,">t.$PREFIX.maximal_mr.names");
#now we process the pb+mega-reads file
foreach $line(@file){
    chomp($line);
    if($line =~ /^>/){
	if(@lines){
	    process_lines(@lines);
	    @lines=();
	}
    }else{
	#print STDERR "Pushed $line\n";
	my @ttt=split(/\s+/,$line);
	push(@lines, \@ttt);
    }
}
#do not forget the last one
process_lines(@lines) if(@lines);

sub process_lines{
    my @args=@_;
    my $o=mummer::Options->new;
    $o->minmatch(10);
    $o->mincluster(60);
    $o->diagfactor(0.2);
    $o->maxgap(200);
    $o->breaklen(120);
    $o->forward();
    my $seq="";
    my $sum_chunk_size=0;
    my $num_chunks=0;
    my $pb_offset=0;
    my $mr_offset=0;
    my $slack=200;

    for(my $i=0;$i<=$#args;$i++){
    
       my ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
       next if($mbgn>$mend);
       next if($bgn>$end);
       #$pb_offset=$bgn-int($mbgn*1.2);
       #$pb_offset=0 if($pb_offset<0);
       #$lpb=($mlen-$mend)*1.2+$end-$pb_offset;
       #$lpb=length($pbseq{$pb})-$pb_offset if($lpb+$pb_offset>length($pbseq{$pb}));
       #$lmr=$mlen;
       #print STDERR "received $bgn,$end,$mbgn,$mend,$mlen,$pb\n";
       if($bgn>$slack){
	$pb_offset=$bgn-$slack-1;
	}else{
	$pb_offset=0;
	}	
	if($mbgn>$slack){
	$mr_offset=$mbgn-$slack-1;
	}else{
	$mr_offset=0;
	}

       if(not defined($readnames{$name})){
	$readnames{$name}=$readnumber;
	print OUTFILE1 ">$readnumber\n$mseq\n";
	print OUTFILE2 "$name\n";
        $readnumber++;
       }
       #print "$bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name\n";
       $lpb=$pb_offset>0?$end-$bgn+2*$slack:$end+$slack;
       $lpb=length($pbseq{$pb})-$pb_offset-1 if($lpb+$pb_offset>length($pbseq{$pb}));
       $lmr=$mr_offset>0?$mend-$mbgn+2*$slack:$mend+$slack;
       $lmr=$mlen-$mr_offset-1 if($lmr+$mr_offset>$mlen);
       #print STDERR "Refine $bgn,$end,$mbgn,$mend,$mlen,$pb\n";
       my $a = mummer::align_sequences(substr($pbseq{$pb},$pb_offset,$lpb), substr($mseq,$mr_offset,$lmr), $o);
       for($j=0;$j<@$a;$j++){

        #print STDERR ">$pb $readnames{$name} ",length($pbseq{$pb})," $mlen\n",$$a[$j]{sA}+$pb_offset," ",$$a[$j]{eA}+$pb_offset," ",$$a[$j]{sB}+$mr_offset," ",$$a[$j]{eB}+$mr_offset," $$a[$j]{Errors} $$a[$j]{SimErrors} $$a[$j]{NonAlphas}\n";
	print ">$pb $readnames{$name} ",length($pbseq{$pb})," $mlen\n",$$a[$j]{sA}+$pb_offset," ",$$a[$j]{eA}+$pb_offset," ",$$a[$j]{sB}+$mr_offset," ",$$a[$j]{eB}+$mr_offset," $$a[$j]{Errors} $$a[$j]{SimErrors} $$a[$j]{NonAlphas}\n";
	print "0\n";
	}
	}
    return();
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

