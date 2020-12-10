#!/usr/bin/env perl
#this code corrects long reads with k-unitigs
use strict;
use warnings;
#first we load the long reads
my $kmer=19;
$kmer=$ARGV[0] if(defined($ARGV[0]));
my $min_density=0.0;
$min_density=$ARGV[1] if(defined($ARGV[1]));
my $DEBUG=0;
my %seq=();
my $seq="";
my $rn="";
my $readname="";
my $seq_local="";
my $line;
my @corrects=();

#now we read the create_mega_reads output and substitute k-unitigs where they match
while($line=<STDIN>){
  chomp($line);
  if(substr($line,0,1) eq ">"){
    if(not($readname eq "")){
      do_corrects(\@corrects) if(length($seq_local)>500);
    }
    ($readname,$seq_local)=split(/\s+/,$line);
    $seq_local=lc($seq_local);
    $readname=substr($readname,1);
    @corrects=();
  }else{
    my @f=split(/\s+/,$line);
#1478.00 1518.00 1479 1499 2 22 2 0.1000 1139641F 41 CAAGATCTGGTCAGAAGAAAAGTTACGTATAGGAATGGATT
    my $lr_start=$f[2]-$f[4]+1;
    my $lr_end=$f[3]+($f[9]-$f[5]);
    my $lr_start_adj=$lr_start;
    $lr_start_adj=1 if($lr_start<1);
    my $lr_end_adj=$lr_end;
    $lr_end_adj=length($seq_local) if($lr_end_adj>length($seq_local));
    my $density=($f[5]-$f[4]+1)/($lr_end_adj-$lr_start_adj+1); #density is calculated over the implied natching region
    my @a=($lr_start,$lr_end,$f[10],$f[8],$density);
    push(@corrects,\@a) if($density>$min_density);
  }
}
$seq{$readname}="";
do_corrects(\@corrects) if(length($seq_local)>500);;

sub do_corrects{
#here we do the actual correction
  my @corrects_sorted=sort {$b->[0] <=> $a->[0]} @{$_[0]};
  my $prevseq="";
  my $prevstart=1000000000;
  my $prevend=1000000000;
  my $prevname="";
  my %bad_k_unitigs=();
  my %overlaps=();
  #first pass we look for bad k-unitigs by examining the overlaps
  #an implied overlap without an actual overlap means that one of the k-unitig alignments is spurious
  #we can check the overlaps on both sides to figure out which one is spurious
  #overlaps are coded as follows
  #number means an actual overlap
  #-1 means an implied overlap without a real overlap
  #0 means no overlap
  print "pass1\n" if($DEBUG);  
  my %overlap_before=();
  my %overlap_after=();
  my $last_overlap=0;
  foreach my $c (@corrects_sorted){
    my $lr_start=$c->[0];
    my $lr_end=$c->[1];
    my $kseq=$c->[2];
    my $kname=$c->[3];
    my $density=$c->[4];
    print "$lr_start $lr_end $kseq $kname $density\n" if($DEBUG);
    
    if($lr_end>$prevend){
      $bad_k_unitigs{$kname}=1;
      print "BAD kunitig non-maximal $kname\n" if($DEBUG);
      next;
    }
    $overlap_before{$kname}=0;
    $overlap_after{$prevname}=0;
    if($prevstart<=$lr_end){#look for an exact overlap of less than 19bp and more than 6bp
      my $overlap=$lr_end-$prevstart+1;
      my $min_overlap=5;
      $min_overlap=$overlap-2 if($overlap<=$min_overlap);
      $min_overlap=1 if($min_overlap<1);
      print "overlap= $overlap, $lr_start $lr_end $prevstart\n" if($DEBUG);
      my $i=0;
      for($i=$kmer-1;$i>$min_overlap;$i--){ #here we test for overlap that is equal or less than the implied overlap and equal or less than k-1
        my $offset=index($kseq,substr($prevseq,0,$i),length($kseq)-$i);
        if($offset>-1){
          $overlap_before{$kname}=$i; 
          $overlap_after{$prevname}=$i;
          $overlaps{"$kname $prevname"}=$i;
          print "overlap found: ",$i,", $lr_start $lr_end $prevstart $offset\n" if($DEBUG);
          print "$kseq\n",' 'x$offset,"$prevseq\n" if($DEBUG);
          last;
        }
      }
      if($i<=$min_overlap){
        if($overlap>=$kmer/2){
          print "bad overlap $kname $prevname\n" if($DEBUG);
          $overlap_before{$kname}=-1;
          $overlap_after{$prevname}=-1;
        }
      }
    }
    $prevname=$kname;
    $prevseq=$kseq;
    $prevstart=$lr_start;
    $prevend=$lr_end;
  }
  #now we go though the list again and invalidate k-unitigs that have invalid overlaps on both sides or invalid and no overlap
  foreach my $c (@corrects_sorted){
    my $kname=$c->[3];
    $overlap_before{$kname}=0 if(not(defined($overlap_before{$kname})));
    $overlap_after{$kname}=0 if(not(defined($overlap_after{$kname})));
    print "$kname $overlap_before{$kname} $overlap_after{$kname}\n" if($DEBUG);
    if(($overlap_before{$kname}==-1 && $overlap_after{$kname}==-1)||($overlap_before{$kname}==0 && $overlap_after{$kname}==-1)||($overlap_before{$kname}==-1 && $overlap_after{$kname}==0)){
    print "BAD kunitig match $kname\n" if($DEBUG);
    $bad_k_unitigs{$kname}=1;
    }
  }
  %overlap_before=();
  %overlap_after=();

  #here in pass2 we do the actual corrections
  $prevseq="";
  $prevstart=1000000000;
  print "PASS2\n" if($DEBUG);
  foreach my $c (@corrects_sorted){
    my $lr_start=$c->[0];
    my $lr_end=$c->[1];
    my $kseq=$c->[2];
    my $kname=$c->[3];
    my $density=$c->[4];
    print "$lr_start $lr_end $kseq $kname $density\n" if($DEBUG);
    next if($bad_k_unitigs{$kname});
    my $insert_len=$lr_end-$lr_start+1; #this is the length of sequence to be replaced on the long read -- could be different from the length of kseq
      if($prevstart<=$lr_end){#look for an exact overlap of less than 19bp and more than 6bp
        my $overlap=$lr_end-$prevstart+1;
        my $min_overlap=5;
        $min_overlap=$overlap-2 if($overlap<=$min_overlap);
        $min_overlap=1 if($min_overlap<1);
        print "overlap= $overlap, $lr_start $lr_end $prevstart\n" if($DEBUG);
        if(not(defined($overlaps{"$kname $prevname"}))){
          my $i=0;
          for(my $i=$kmer-1;$i>$min_overlap;$i--){ #here we test for overlap that is equal or less than the implied overlap and equal or less than k-1
            my $offset=index($kseq,substr($prevseq,0,$i),length($kseq)-$i);
            if($offset>-1){
              $insert_len+=($i-$overlap);
              print "overlap found: ",$i,", $lr_start $lr_end $prevstart $offset\n" if($DEBUG);
              print "$kseq\n",' 'x$offset,"$prevseq\n" if($DEBUG);
              last;
            }
          }
          if($i<=$min_overlap){#no overlap -- adjust the interval
            $insert_len-=$overlap;
          }
        }else{
          $insert_len+=($overlaps{"$kname $prevname"}-$overlap);
          print "overlap previously found: ",$overlaps{"$kname $prevname"},", $lr_start $lr_end $prevstart\n" if($DEBUG);
        }
      }
    if($lr_start<1){#if k-unitig extends the beginning of the read
      $seq_local=$kseq.substr($seq_local,$lr_end);
    }elsif($lr_end>length($seq_local)){#if k-unitig extends the end of the read
      $seq_local=substr($seq_local,0,$lr_start).$kseq;
    }else{
      print "inserting $lr_start $insert_len $lr_end\n" if($DEBUG);
      print "before  : ",substr($seq_local,$lr_start-1,length($kseq)+20)," $insert_len ",length($kseq),"\n" if($DEBUG);
      substr($seq_local,$lr_start-1,$insert_len)=$kseq;
      print "after   : ",substr($seq_local,$lr_start-1,length($kseq)+20),"\n" if($DEBUG);
      print "inserted: $kseq\n" if($DEBUG);
    }
    $prevname=$kname;
    $prevseq=substr($seq_local,$lr_start-1,$kmer);
    $prevstart=$lr_start;
  }
  print ">$readname\n$seq_local\n";
}
