#!/usr/bin/env perl
#this code splits reads to join into batches according to join sequences
my $qrys_file=$ARGV[0];
my $prefix=$ARGV[1];

#first we read and assign joining sequence file names to batches
for($i=2;$i<=$#ARGV;$i++){
open(FILE,$ARGV[$i]);
while($line=<FILE>){
if($line=~/^>/){
  chomp($line);
  my @f=split(/\//,$line);
  $batch{substr($f[0],1)}=$i-1; #1-based numbers
}
}
close(FILE);
}
#then we read the file corresponding read names to joining sequence names -- there may be more than one joining sequence per read name!
open(FILE,$qrys_file);
while($line=<FILE>){
chomp($line);
my @f=split(/\s+/,$line);
$joining_name{$f[1]}.="$f[2]_$f[3] ";
}

#now we read the file of reads to join from stdin and for each batch of joining sequences output the corresponding batch of reads.
#ASSUMES read sequences are ONE LINE
for($i=2;$i<=$#ARGV;$i++){
local *OUTFILE;
open(OUTFILE,">$prefix.".($i-1).".fa");
push(@fh,*OUTFILE);
}
while($line=<STDIN>){
if($line=~/^>/){
  chomp($line);
  my @f=split(/\s+/,$line);
  my $readname_j=substr($f[0],1);
  my @ff=split(/\./,$readname_j);
  my $readname=$ff[0];
  #print "DEBUG processing read $readname_j $readname\n";
  $line=<STDIN>;
  chomp($line);
  my $readseq=$line;
  if(defined($joining_name{$readname})){
    my @ff=split(/\s+/,$joining_name{$readname});
    #print "DEBUG ",join(" ",@ff),"\n";
    my %used_batches=();
    foreach my $j(@ff){
      if(defined($batch{$j}) && not(defined($used_batches{$batch{$j}}))){
        #print "DEBUG output $readname into batch $batch{$j}\n";
        print {$fh[$batch{$j}-1]} ">$readname_j\n$readseq\n";
        $used_batches{$batch{$j}}=1;
        }
      }
  }
}
}
