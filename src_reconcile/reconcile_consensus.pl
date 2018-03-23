#!/usr/bin/env perl
########################################
##Copyright Johns Hopkins University 2018#
########################################
##this code refines the alignments of mega reads to pacbio reads using nucmer
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
while($line=<STDIN>){
  print $line if($DEBUG);
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  if(not($f[-2] eq $last_contig)){
    $last_offset=0;
  }
  next if($f[1]<=$last_offset);
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
  my $len_ref=length($rseq{$f[-2]});
  my $adj_beg=$len_ref-($f[11]-$f[0]);
  my $adj_end=$len_ref-($f[11]-$f[1]);
  print "$f[0] $f[1] $adj_beg $adj_end $len_ref\n" if($DEBUG);
  $rseq{$f[-2]}=substr($rseq{$f[-2]},0,$adj_beg-1).$subseq.substr($rseq{$f[-2]},$adj_end);
  $last_contig=$f[-2];
  print length($rseq{$f[-2]}),"\n" if($DEBUG);
  $last_offset=$f[1];
}

foreach $c(keys %rseq){
  print ">$c\n$rseq{$c}\n";
}
