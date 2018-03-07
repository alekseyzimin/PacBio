#!/usr/bin/env perl
########################################
##Copyright Johns Hopkins University 2018#
########################################
##this code refines the alignments of mega reads to pacbio reads using nucmer
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use mummer;
#this code replaces consensus of one assembly with the other based on the nucmer matches -- consensus of the reference is replaced with the query
#
my $DEBUG=0;
my $ref_contigs=$ARGV[0];
my $qry_contigs=$ARGV[1];

my $ctg="",$seq="";
my %rseq =(),%qseq=();
open(FILE,$ref_contigs);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $rseq{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $rseq{$ctg}=$seq;
}

$ctg="",$seq="";
open(FILE,$qry_contigs);
while($line=<FILE>){
  chomp($line);   
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $qseq{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $qseq{$ctg}=$seq;
}



#now read in the coords file. We need this file to figure out where the filtered alignments of the contigs are, but we cannot
#use the alignment coordinates directly due to indels. Therefore we have to re-align the polishing contigs locally
my $subseq;
my $offset=0;
my $last_contig="";
my $padding=1000;
while($line=<STDIN>){
  print $line if($DEBUG);
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  if(not($f[-2] eq $last_contig)){
    $offset=0;
    $last_offset=0;
  }

  next if(not(defined($rseq{$f[-2]})));
  next if(not(defined($qseq{$f[-1]})));
  $subseq="";
  if($f[3]<$f[4]){
    $subseq=substr($qseq{$f[-1]},$f[3]-1,$f[4]-$f[3]+1);
  }else{
    $subseq=substr($qseq{$f[-1]},$f[4]-1,$f[3]-$f[4]+1);
    $subseq=reverse($subseq);
    $subseq=~tr/ACGTNacgtn/TGCAntgcan/;
  }
  my $adj_beg=$f[0]+$offset;
  my $adj_end=$f[1]+$offset;
#here we cut-out the part of the reference between the new coordinates padded by 1000 bp on each side and re-align
  my $ref_sub_beg=$adj_beg-$padding;
  $ref_sub_beg=0 if($ref_sub_beg<0);
  print "ref_sub starts at $ref_sub_beg $adj_beg $adj_end\n" if($DEBUG);
  my $ref_to_align=substr($rseq{$f[-2]},$ref_sub_beg,($adj_end-$adj_beg)+$padding);
  #calling nucmer, minimum match cluster of 200 because we are not interested in short alignments
  my $o=mummer::Options->new;
  $o->minmatch(31);
  $o->mincluster(100);
  $o->forward();
  my $a = mummer::align_sequences($ref_to_align,$subseq,$o);
  print "found matches:",scalar(@$a),"\n" if($DEBUG);
  #pick the longest match
  my $lindex=-1;
  my $matchlen=0;
  for($j=0;$j<@$a;$j++){
    if($$a[$j]{eB}-$$a[$j]{sB}>$matchlen){
      $matchlen=$$a[$j]{eB}-$$a[$j]{sB};
      $lindex=$j;
    }
    print $$a[$j]{sA}+$ref_sub_beg," ",$$a[$j]{eA}+$ref_sub_beg," ",$$a[$j]{sB}," ",$$a[$j]{eB},"\n" if($DEBUG);
  }
  print "picked match $lindex\n" if($DEBUG);
  next if($lindex==-1);#wrong orientation
  $subseq_adj=substr($subseq,$$a[$lindex]{sB}-1,$$a[$lindex]{eB}-$$a[$lindex]{sB}+1);
  $rseq{$f[-2]}=substr($rseq{$f[-2]},0,$$a[$lindex]{sA}+$ref_sub_beg-1).$subseq_adj.substr($rseq{$f[-2]},$$a[$lindex]{eA}+$ref_sub_beg);
  $offset+=(length($subseq)-($adj_end-$adj_beg)-1);
  $last_contig=$f[-2];
  print length($rseq{$f[-2]}),"\n" if($DEBUG);
}

foreach $c(keys %rseq){
  print ">$c\n$rseq{$c}\n";
}
