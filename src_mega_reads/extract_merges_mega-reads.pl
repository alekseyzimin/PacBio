#!/usr/bin/env perl
#
#this code extracts possible contig merges from a nucmer alignment: the reference sequences are merged with query sequence intput is a delta file 
#ASSUMES show-coords -q output, that is sorted by query coord!!!!!!

#open(FILE,"delta-filter -q -i 98 $ARGV[0] | show-coords -lcHq -L 1000 /dev/stdin |"); 
my $slack=400;
my $maxgap=100000;
my $mingap=-50;
$slack = $ARGV[2] if(not($ARGV[2] eq ""));
$mingap = $ARGV[3] if(not($ARGV[3] eq ""));
$maxgap = $ARGV[4] if(not($ARGV[4] eq ""));

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

open(FILE,$ARGV[1]);#file valid merges 
while($line=<FILE>){
  chomp($line);
  $valid_merges{$line}=1;
}

#first we read in all the matches into an array
my $prevline="";
while($line=<STDIN>){
    chomp($line);
    $line=~s/^\s+//;
    my @f=split(/\s+/,$line);
    push(@lines,$line) if($f[0]<=$slack || $f[1]>=$f[11]-$slack);#only see potentially useful lines
}

#now we go through the array and collect the possible merges
#
#
for($i=0;$i<=$#lines;$i++){
  @f1=split(/\s+/,$lines[$i]);
  ($r1,$junk)=split(/\./,$f1[-2]);
  for($j=$i+1;$j<=$#lines;$j++){ 
    @f2=split(/\s+/,$lines[$j]);
    #if switched to the next reference contig, stop
    if(not($f2[-1] eq $f1[-1])){$j=$#lines;next;}
    ($r2,$junk)=split(/\./,$f2[-2]);
    #next if matches to the same contig or chunks of different contigs
    next if($f1[-2] eq $f2[-2] || not($r1 eq $r2));
    #next if merge not valid
    next if(not(defined($valid_merges{"$f1[-2] $f2[-2]"})) && not(defined($valid_merges{"$f2[-2] $f1[-2]"})));
    #print STDERR "testing\n$lines[$i]\n$lines[$j]\n";
    #now testing merge
    if($f1[3]<$f1[4]){
      if($f2[3]<$f2[4]){
      #forward merge ---->     ------>
        $gap=$f2[3]-$f1[4];
        my $trim_e=$f1[11]-$f1[1];
        my $trim_b=$f2[0]-1;
        if($trim_e<=$slack && $trim_b<=$slack && $gap<$maxgap && $gap>$mingap && $valid_merges{"$f1[-2] $f2[-2]"}){
          #$valid_merges{"$f1[-2] $f2[-2]"}=0;
          print "$f1[-2] $trim_e F $f2[-2] $trim_b F $gap ";
          die("Query sequence $f1[-1] is not found") if(not(defined($qseq{$f1[-1]})));
          print substr($qseq{$f1[-1]},$f1[4],$gap) if($gap>0);
          print "\n";
          #last;
        }
      }
    }else{
      if($f2[3]>$f2[4]){
      #reverse merge  <---------   <----------
        $gap=$f2[4]-$f1[3];
        my $trim_e=$f1[0]-1;
        my $trim_b=$f2[11]-$f2[1];
        if($trim_e<$slack && $trim_b<$slack && $gap<$maxgap  && $gap>$mingap &&  $valid_merges{"$f2[-2] $f1[-2]"}){
          #$valid_merges{"$f2[-2] $f1[-2]"}=0;
          print "$f1[-2] $trim_e R $f2[-2] $trim_b R $gap ";
          die("Query sequence $f1[-1] is not found") if(not(defined($qseq{$f1[-1]})));
          print substr($qseq{$f1[-1]},$f1[3],$gap) if($gap>0);
          print "\n";
          #last;
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
