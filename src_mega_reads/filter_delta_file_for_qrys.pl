#!/usr/bin/env perl
#this code filters the matches in delta file for alignments contained in qrys.txt
my $qrys_file=$ARGV[0];
open(FILE,$qrys_file);
while($line=<FILE>){
  chomp($line);
  @f=split(/\s+/,$line); 
  $h{"$f[1] $f[2]_$f[3]"}=1;
  }
#now parse delta file on stdin -- first two lines are the header
$line=<STDIN>;
print $line;
$line=<STDIN>;
print $line;
$output=0;
while($line=<STDIN>){
  if($line =~ /^>/){
    chomp($line);
    @f1=split(/\s+/,substr($line,1));
    @f2=split(/\//,$f1[1]); 
    @f3=split(/\./,$f1[0]); 
    if(defined($h{"$f3[0] $f2[0]"})){
      $output=1; 
      $hline=$line;
      $houtput=1;
    }else{
      $output=0;
    }
  }elsif($output){
    chomp($line);
    @f4=split(/\s+/,$line);
    if(scalar(@f4)>1){
      if($houtput){
        print "$hline\n";
        $houtput=0;
      } 
      print "$line\n0\n";
    }
  }
} 
