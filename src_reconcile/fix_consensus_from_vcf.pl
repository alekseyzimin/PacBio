#!/usr/bin/env perl
#this code fixes errors in the consensus called in vcf file by freebayes
#first we read the sequences into memory
#
my $min_alt_obs=2;
my $ref_contigs=$ARGV[0];
my $ctg="",$seq="";
my %rseq =();
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

#we assume the vcf file is sorted by contig+coordinate; read in all vcf entries for the first contig
$ctg="";
my @fixes=();
my @originals=();
my $offsets=();
while($line=<STDIN>){
next if($line =~ /^\#/);
my @f=split(/\s+/,$line);
if(not($f[0] eq $ctg)){
  if($#fixes>=0){
  #proceed with fixing
    for(my $i=$#fixes;$i>-1;$i--){#going in reverse order to avoid shifting sequence due to indels
      die("sequence $ctg not found in the input fasta file") if(not(defined($rseq{$ctg})));
      #first we check if the sequence at given offset matches the original variant
      my $original_seq=substr($rseq{$ctg},$offsets[$i]-1,length($originals[$i]));
      #print "$i $ctg $offsets[$i] $originals[$i] $original_seq $fixes[$i]\n";
      die("sequence does not match the original $original_seq $originals[$i]") if(not($originals[$i] eq $original_seq));
      #then substitute
      $rseq{$ctg}=substr($rseq{$ctg},0,$offsets[$i]-1).$fixes[$i].substr($rseq{$ctg},$offsets[$i]+length($originals[$i])-1);
      }
  }
  @fixes=();
  @originals=();
  @offsets=();
  $ctg=$f[0];
}
my @ff=split(/:/,$f[9]);
#print "$f[9]\n";
if($ff[3]==0 && $ff[5]>=$min_alt_obs){
  push(@fixes,$f[4]);
  push(@originals,$f[3]);
  push(@offsets,$f[1]);
}
}

if($#fixes>=0){
#proceed with fixing
  for(my $i=$#fixes;$i>-1;$i--){#going in reverse order to avoid shifting sequence due to indels
#first we check if the sequence at given offset matches the original variant
    die("sequence $ctg not found in the input fasta file") if(not(defined($rseq{$ctg})));
    my $original_seq=substr($rseq{$ctg},$offsets[$i]-1,length($originals[$i]));
    #print "$offsets[$i] $originals[$i] $original_seq $fixes[$i]\n";
    die("sequence does not match the original $original_seq $originals[$i]") if(not($originals[$i] eq $original_seq));
#then substitute
    $rseq{$ctg}=substr($rseq{$ctg},0,$offsets[$i]-1).$fixes[$i].substr($rseq{$ctg},$offsets[$i]+length($originals[$i])-1);
  }
}

foreach my $c(keys %rseq){
  print ">$c\n$rseq{$c}\n";
}

