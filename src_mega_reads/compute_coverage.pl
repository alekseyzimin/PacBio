#!/usr/bin/env perl
#
#This code takes as input file
#read ctg start
#read ctg end
#sorted by the third column and returns the coverage for ctg
#
#
my $cctg="";
my %reads=();
my $coverage=0;
while($line=<STDIN>){
chomp($line);
($read,$ctg,$pos)=split(/\s+/,$line);
if($ctg eq $cctg){
if(defined($reads{$read})){
$coverage--;
undef($reads{$read});
}else{
$coverage++;
$reads{$read}=1;
}
}else{
$coverage=1;
%reads=();
$reads{$read}=1;
$cctg=$ctg;
}
print "$line $coverage\n"
}
