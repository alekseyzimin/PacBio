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
my $slack=1000;
my $maxgap=10000;
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
for($i=0;$i<$#lines;$i++){
  @f1=split(/\s+/,$lines[$i]);
  $j=$i+1;
#print "DEBUG $i $j\n$lines[$i]\n$lines[$j]\n\n"; 
    @f2=split(/\s+/,$lines[$j]);
    next if($f1[-2] eq $f2[-2]);
    next if(not($f1[-1] eq $f2[-1]));
    my $oh1=0;
    my $oh2=0;
    my $success=0;
    my $gstart=1;
#print "DEBUG considering $i $j\n$lines[$i]\n$lines[$j]\n\n";
    if($f1[3]<$f1[4]){
      $gstart=$f1[4];
      if($f2[3]<$f2[4]){
#forward forward merge ---->     ------>
        $gap=$f2[3]-$f1[4]+1;
        $oh1=$f1[11]-$f1[1];
        $oh2=$f2[0]-1;
        if($oh1<$slack && $oh2<$slack && $gap<$maxgap && $gap>$mingap){
          $success=1;
          if($f1[-2] lt $f2[-2]){
            print "$f1[-2] $oh1 F $f2[-2] $oh2 F $gap ",;
          }else{
            print "$f2[-2] $oh2 R $f1[-2] $oh1 R $gap ";
          }
        }
      }else{
#forward reverse merge ---->     <------
        $gap=$f2[4]-$f1[4]+1;
        $oh1=$f1[11]-$f1[1];
        $oh2=$f2[11]-$f2[1];
        if($oh1<$slack && $oh2<$slack && $gap<$maxgap && $gap>$mingap){
          $success=1;
          if($f1[-2] lt $f2[-2]){
            print "$f1[-2] $oh1 F $f2[-2] $oh2 R $gap ";
          }else{
            print "$f2[-2] $oh2 F $f1[-2] $oh1 R $gap ";
          }
        }
      }
    }else{
      $gstart=$f1[3];
      if($f2[3]<$f2[4]){
#reverse forward merge <-----     ------>
        $gap=$f2[3]-$f1[3]+1;
        $oh1=$f1[0]-1;
        $oh2=$f2[0]-1;
        if($oh1<$slack && $oh2<$slack && $gap<$maxgap && $gap>$mingap){
          $success=1;
          if($f1[-2] lt $f2[-2]){
            print "$f1[-2] $oh1 R $f2[-2] $oh2 F $gap ";
          }else{
            print "$f2[-2] $oh2 R $f1[-2] $oh1 F $gap ";
          }
        }
      }else{
#reverse reverse merge <-----     <------
        $gap=$f2[4]-$f1[3]+1;
        $oh1=$f1[0]-1;
        $oh2=$f2[11]-$f2[1];
        if($oh1<$slack && $oh2<$slack && $gap<$maxgap && $gap>$mingap){
          $success=1;
          if($f1[-2] lt $f2[-2]){
            print "$f1[-2] $oh1 R $f2[-2] $oh1 R $gap ";
          }else{
            print "$f2[-2] $oh2 F $f1[-2] $oh2 F $gap ";
          }
        } 
      }
    }
    if($success){
      if($gap>0){
        $gstart=1 if($gstart<1);
        if($f1[-2] lt $f2[-2]){
          print lc(substr($qseq{$f1[-1]},$gstart-1,$gap));
        }else{
          print reverse_complement(lc(substr($qseq{$f1[-1]},$gstart-1,$gap)));
        }
      }else{
        print "n";
      }
      print " $sum_overhangs\n";
      $success=0;
    }
  #}#$j loop
}#$i loop

sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}
