#!/usr/bin/env perl
#this code translates into reduced read names and positions
#
my $reduceFile =$ARGV[0];
open(FILE,$reduceFile);
while($line=<FILE>){
chomp($line);
my ($containee, $container, $ori, $offset)=split(/\s+/,$line);
$reduced{$containee}="$container $ori $offset";
}

while($line=<STDIN>){
chomp($line);
my ($rname,$srname,$offset,$ori,$status)=split(/\s+/,$line);
if(defined($reduced{$srname})){
  my($container, $cori, $coffset)=split(/\s+/,$reduced{$srname});
  if($cori eq "F"){
    $offset+=$coffset;
  }else{
    $ori=~tr/FR/RF/;
    $offset=$coffset-$offset;
  }
  print "$rname $container $offset $ori\n";
}else{
  print "$rname $srname $offset $ori\n";
}
}
