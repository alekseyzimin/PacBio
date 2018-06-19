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
my $slack=500;
my $maxgap=100000;
my $mingap=-10;
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
  ($r1,$junk)=split(/\./,$f1[-2]);
  for($j=$i+1;$j<=$#lines;$j++){ 
#print "DEBUG $i $j\n$lines[$i]\n$lines[$j]\n\n"; 
    @f2=split(/\s+/,$lines[$j]);
    ($r2,$junk)=split(/\./,$f2[-2]);
    next if($f1[-2] eq $f2[-2] || not($f1[-1] eq $f2[-1]) || not($r1 eq $r2));
    if($f1[3]<$f1[4]){
      if($f2[3]<$f2[4]){
      #forward join ---->     ------>
        $gap=$f2[3]-$f1[4];
        my $trim_e=$f1[11]-$f1[1];
        my $trim_b=$f2[0]-1;
        if($trim_e<=$slack && $trim_b<=$slack && $gap<$maxgap && $gap>$mingap){
          print "$f1[-2] $trim_e F $f2[-2] $trim_b F $gap ";
          print substr($qseq{$f1[-1]},$f1[4],$gap) if($gap>0);
          print "\n";
          last;
        }
      }
    }else{
      if($f2[3]>$f2[4]){
      #reverse join   <---------   <----------
        $gap=$f2[4]-$f1[3];
        my $trim_e=$f1[0]-1;
        my $trim_b=$f2[11]-$f2[1];
        if($trim_e<$slack && $trim_b<$slack && $gap<$maxgap  && $gap>$mingap){
          print "$f1[-2] $trim_e R $f2[-2] $trim_b R $gap ";
          print substr($qseq{$f1[-1]},$f1[3],$gap) if($gap>0);
          print "\n";
          last;
        } 
      }
    }
  }
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}
