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
	if(@lines && $#lines<100){
	    process_lines(@lines);#no more than 100 chunks per PB read
	    @lines=();
	}
	($rn,$junk)=split(/\s+/,substr($line,1));
    }else{
	my @ttt=split(/\s+/,$line);
	push(@lines, \@ttt);
    }
}
#do not forget the last one
process_lines(@lines) if(@lines);

sub process_lines{
    my @args=@_;
    my $o=mummer::Options->new;
    $o->minmatch(15);
    $o->mincluster(100);
    $o->diagfactor(0.2);
    $o->maxgap(200);
    $o->breaklen(100);
    $o->forward();
    my $seq="";
    my $sum_chunk_size=0;
    my $num_chunks=0;

    for(my $i=0;$i<=$#args;$i++){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
        $sum_chunk_size+=$mlen;
        $num_chunks++;
        }

    return() if($sum_chunk_size/$num_chunks<200);#average chunk size must be >200bp

    for(my $i=0;$i<=$#args;$i++){
       my ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
       if(not defined($readnames{$name})){
	$readnames{$name}=$readnumber;
	print OUTFILE1 ">$readnumber\n$mseq\n";
	print OUTFILE2 "$name\n";
        $readnumber++;
       }
       #print "$bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name\n";
       my $a = mummer::align_sequences($pbseq{$pb}, $mseq, $o);
       for($j=0;$j<@$a;$j++){
	print ">$pb $readnames{$name} ",length($pbseq{$pb})," $mlen\n";
	print "$$a[$j]{sA} $$a[$j]{eA} $$a[$j]{sB} $$a[$j]{eB} $$a[$j]{Errors} $$a[$j]{SimErrors} $$a[$j]{NonAlphas}\n";
#foreach my $d($$a[$j]{delta}){
#		print $d,"\n";
#	}
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

