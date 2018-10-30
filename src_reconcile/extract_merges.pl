#!/usr/bin/env perl
#
#this code extracts possible contig merges from a nucmer alignment: the reference sequences are merged with query sequence intput is a delta file 
#ASSUMES show-coords -q output, that is sorted by query coord!!!!!!

#open(FILE,"delta-filter -q -i 98 $ARGV[0] | show-coords -lcHq -L 1000 /dev/stdin |"); 

open(FILE,$ARGV[0]);#file with query contigs
while($line=<FILE>){
    chomp($line);
    if($line =~ /^>/){
	@f=split(/\s+/,$line);
	$qn=substr($f[0],1);
    }else{
	$qseq{$qn}.=$line;
    }
}
#first we read in all the matches into an array
my $prevline="";
my $slack=5;
my $maxgap=100000;
my $mingap=-10000;
while($line=<STDIN>){
    chomp($line);
    $line=~s/^\s+//;
    my @f=split(/\s+/,$line);
    push(@lines,$line) if($f[0]<=$slack || $f[1]>=$f[11]-$slack);#only see useful lines
}

#now we go through the array and collect the possible merges
#
#
for($i=0;$i<=$#lines;$i++){
@f1=split(/\s+/,$lines[$i]);
for($j=$i+1;$j<=$#lines;$j++){ 
#print "DEBUG $i $j\n$lines[$i]\n$lines[$j]\n\n"; 
    @f2=split(/\s+/,$lines[$j]);
    next if($f1[-2] eq $f2[-2]);
    next if(not($f1[-1] eq $f2[-1]));
	if($f1[3]<$f1[4]){
	    if($f2[3]<$f2[4]){
		$gap=$f2[3]-$f1[4];
		if($f1[1]>$f1[11]-$slack && $f2[0]<$slack && $gap<$maxgap && $gap>$mingap){
		    print "$f1[-2] F $f2[-2] F $gap ";
		    print substr($qseq{$f1[-1]},$f1[4],$gap) if($gap>0);
		    print "\n";
                    last;
		}
	    }else{
		$gap=$f2[4]-$f1[4];
		if($f1[1]>$f1[11]-$slack && $f2[1]>$f2[11]-$slack && $gap<$maxgap  && $gap>$mingap){
		    print "$f1[-2] F $f2[-2] R $gap ";
		    print substr($qseq{$f1[-1]},$f1[4],$gap) if($gap>0);
		    print "\n";
                    last;
		}
	    }
	}else{
	    if($f2[3]<$f2[4]){
		$gap=$f2[3]-$f1[3];
		if($f1[0]<$slack && $f2[0]<$slack && $gap<$maxgap  && $gap>$mingap){
		    print "$f1[-2] R $f2[-2] F $gap ";
		    print substr($qseq{$f1[-1]},$f1[3],$gap) if($gap>0);
		    print "\n";
                    last;
		}
	    }else{
		$gap=$f2[4]-$f1[3];
		if($f1[0]<$slack && $f2[1]>$f2[11]-$slack && $gap<$maxgap  && $gap>$mingap){
		    print "$f1[-2] R $f2[-2] R $gap ";
		    print substr($qseq{$f1[-1]},$f1[3],$gap) if($gap>0);
		    print "\n";
                    last;
		} 
	    }
	}
        #last if($gap>$maxgap);
}
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}
