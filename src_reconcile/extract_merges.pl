#!/usr/bin/env perl
#
#this code extracts possible contig merges from a nucmer alignment: the reference sequences are merged with query sequence intput is a delta file 
#ASSUMES show-coords -q output, that is sorted by query coord!!!!!!

my $output_prefix="patches";
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
my $slack=1000;
my $maxgap=200000;
my $mingap=-10000;
my %oh1=(),%oh2=(),%gap=(),%gapseq=();
while($line=<STDIN>){
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  push(@lines,$line) if($f[0]<=$slack || $f[1]>=$f[11]-$slack);#only see useful lines
}

#now we go through the array and collect the possible merges
#
#
for($i=0;$i<$#lines;$i++){
  @f1=split(/\s+/,$lines[$i]);
  $j=$i+1;
#print "DEBUG $i $j\n$lines[$i]\n$lines[$j]\n\n"; 
    @f2=split(/\s+/,$lines[$j]);
    next if($f1[-2] eq $f2[-2]);
    next if(not($f1[-1] eq $f2[-1]));
    my $oh1=0;
    my $oh2=0;
    my $gstart=1;
    my $dir1="F";
    my $dir2="F";
#print "DEBUG considering $i $j\n$lines[$i]\n$lines[$j]\n\n";
    if($f1[3]<$f1[4]){
      $gstart=$f1[4];
      if($f2[3]<$f2[4]){
#forward forward merge ---->     ------>
        $gap=$f2[3]-$f1[4]+1;
        $oh1=$f1[11]-$f1[1];
        $oh2=$f2[0]-1;
        $dir1="F";
        $dir2="F";
      }else{
#forward reverse merge ---->     <------
        $gap=$f2[4]-$f1[4]+1;
        $oh1=$f1[11]-$f1[1];
        $oh2=$f2[11]-$f2[1];
        $dir1="F";
        $dir2="R";
      }
    }else{
      $gstart=$f1[3];
      if($f2[3]<$f2[4]){
#reverse forward merge <-----     ------>
        $gap=$f2[3]-$f1[3]+1;
        $oh1=$f1[0]-1;
        $oh2=$f2[0]-1;
        $dir1="R";
        $dir2="F";
      }else{
#reverse reverse merge <-----     <------
        $gap=$f2[4]-$f1[3]+1;
        $oh1=$f1[0]-1;
        $oh2=$f2[11]-$f2[1];
        $dir1="R";
        $dir2="R";
        } 
    }
    if($oh1<$slack && $oh2<$slack && $gap<$maxgap && $gap>$mingap){
        $gstart=1 if($gstart<1);
        if($f1[-2] lt $f2[-2]){
          $joinline="$f1[-2]:$dir1:$f2[-2]:$dir2";
          #print "DEBUG $joinline $oh1 $oh2 $gap\n";
          if(not(defined($oh1{$joinline})) || $oh1{$joinline} + $oh2{$joinline} > $oh1+$oh2){
            $gseq{$joinline}="n";
            $gseq{$joinline}=lc(substr($qseq{$f1[-1]},$gstart-1,$gap)) if($gap>0);
            $oh1{$joinline}=$oh1;
            $oh2{$joinline}=$oh2;
            $gap{$joinline}=$gap;
          }
          $paircount{"$f1[-2] $f2[-2]"}++;
        }else{
          $dir1= $dir1 eq "F" ? "R" : "F";
          $dir2= $dir2 eq "F" ? "R" : "F";
          $joinline="$f2[-2]:$dir2:$f1[-2]:$dir1";
          #print "DEBUG $joinline $oh1 $oh2 $gap\n";
          if(not(defined($oh1{$joinline})) || $oh1{$joinline} + $oh2{$joinline} > $oh1+$oh2){
            $gseq{$joinline}="n";
            $gseq{$joinline}=reverse_complement(lc(substr($qseq{$f1[-1]},$gstart-1,$gap))) if($gap>0);
            $oh1{$joinline}=$oh2;
            $oh2{$joinline}=$oh1;
            $gap{$joinline}=$gap;
          }
          $paircount{"$f2[-2] $f1[-2]"}++;
        }
      $joincount{$joinline}++;
      $rseq{$joinline}.=$qseq{$f1[-1]}." ";
    }
}#$i loop

#now we have all information in the hashes, let's output the bundles
#the longest read is the seq, the rest used to polish
#we output the reads and run the consensus for each patch separately
foreach my $k (keys %rseq){
  my @seq=split(/\s+/,$rseq{$k});
  my $max_len=0;
  my $max_i=0;
  for($i=0;$i<=$#seq;$i++){
    if(length($seq[$i])>$max_len){
      $max_i=$i;
      $max_len=length($seq[$i]);
    }
  }
  $seqs{$seq[$max_i]}=1;
  for($i=0;$i<=$#seq;$i++){
    $seqs{$seq[$max_i]}.=" ".$seq[$i] unless($i==$max_i);
  }
}
#flatten the seqs and run consensus
#do_consensus.sh must exist!

if(-e "do_consensus.sh"){
  open(RAW,">patches.raw.fa");
  my $index=0;
  foreach my $seq(keys %seqs){
    $index++;
    my @seq=split(/\s+/,$seqs{$seq});
    if($#seq==0){#no polishing seq -- put into raw
      print RAW ">$index\n$seq\n";
    }else{
      open(REF,">patches.ref.fa");
      open(READS,">patches.reads.fa");
      print REF ">$index\n$seq\n";
      #uniq the consensus reads
      my %uniqh=();
      for(my $i=1;$i<=$#seq;$i++){
        $uniqh{$seq[$i]}=1;
      }
      my $index_r=0;
      foreach my $k(keys %uniqh){
        print READS ">r$index_r\n$k\n";
        $index_r++;
      }
      close(REF);
      close(READS);
#run consensus
      system("./do_consensus.sh");
    }
  }
}

#output the links
foreach $k (keys %joincount){
  my @f=split(/:/,$k);
  if($paircount{"$f[0] $f[2]"} == $joincount{$k} || $joincount{$k}>1){
    print "$f[0] $oh1{$k} $f[1] $f[2] $oh2{$k} $f[3] $gap{$k} $gseq{$k}\n";
  }
}

sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}
