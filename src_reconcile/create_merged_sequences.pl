#!/usr/bin/env perl
#first we read in sequences of contigs to merge
my $contigs_to_merge=$ARGV[0];
my $ctg="",$seq="";
my %seq =(),%output=();
open(FILE,$contigs_to_merge);
while($line=<FILE>){
    chomp($line);
    if($line=~/^\>/){
	my @f=split(/\s+/,$line);
	if(not($seq eq "")){
	    $seq{$ctg}=$seq;
	}
	$ctg=substr($f[0],1);
	$seq="";
    }else{
	$seq.=$line;
    }
}
if(not($seq eq "")){
    $seq{$ctg}=$seq;
}

#then we read in the merging sequences; we record the sequence that follows each contig
my $merges=$ARGV[1];
open(FILE,$merges);
while($line=<FILE>){
    #print STDERR $line;
    chomp($line);
    my ($c1,$d1,$c2,$d2,$g,$s)=split(/\s+/,$line);
    $gseq{"$c1$d1$c2$d2"}=$s;
    $d1 =~ tr/FR/RF/;
    $d2 =~ tr/FR/RF/;
    $gseq{"$c2$d2$c1$d1"}=reverse_complement($s);
} 


#now read in the merges file
while($line=<STDIN>){
    chomp($line);
    my @f=split(/\s+/,$line);
    print ">",join("_",@f),"\n";
    #print STDERR ">",join("_",@f),"\n";
#output first contig
    if($f[1] eq "R"){
	print reverse_complement($seq{$f[0]});
    }else{
	print $seq{$f[0]};
    }
    $output{$f[0]}=1;
#now the same for the rest we first output the previous gap (or trim) and then the contig
    for($i=3;$i<$#f;$i+=3){
    #print STDERR "$i\n";
	if($f[$i-1]>0){
            die("gap $f[$i-3]$f[$i-2]$f[$i]$f[$i+1] not found") if(not(defined($gseq{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"})));
            die("sequence $f[$i] not found") if(not(defined($seq{$f[$i]})));
            print $gseq{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"};
            #print STDERR "gap $f[$i-3]$f[$i-2]$f[$i]$f[$i+1] ",$gseq{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"},"\n";
	    if($f[$i+1] eq "R"){
		print reverse_complement($seq{$f[$i]});
	    }else{
		print $seq{$f[$i]};
	    }
	}else{#negative gap
	    if($f[$i+1] eq "R"){
		print substr(reverse_complement($seq{$f[$i]}),-$f[$i+2]);
	    }else{
		print substr($seq{$f[$i]},-$f[$i+2]);
	    }
	}
	$output{$f[$i]}=1;
    }
    print "\n";
}

foreach $c(keys %seq){
    print ">$c\n$seq{$c}\n" if(not(defined($output{$c})));
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

