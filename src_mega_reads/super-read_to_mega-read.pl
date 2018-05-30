#!/usr/bin/env perl
#
#this code, given the output of the super-mega-reads, eliminate contained mega-reads
my $readPlacementFile=$ARGV[0];
my $megaReadsFile=$ARGV[1];
my $superReadNamesFile=$ARGV[2];
my @mr_sizes;
my @mr_names;
my %groups;

open(FILE,$megaReadsFile);
while($line=<FILE>){
if(substr($line,0,1) eq ">"){
chomp($line);
push(@mr_names,substr($line,1));
}else{
push(@mr_sizes,length($line));
}
}

open(FILE,$superReadNamesFile);
while($line=<FILE>){
chomp($line);
push(@sr_names,$line);
}

open(FILE,$readPlacementFile);
while($line=<FILE>){
chomp($line);
@f=split(/\s+/,$line);
my $mrn=int(substr($f[0],2)/2);
print "$mr_names[$mrn] $f[3] $sr_names[$f[1]]\n";
}
