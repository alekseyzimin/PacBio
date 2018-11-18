#!/usr/bin/env perl
#this code translates into reduced read names and positions
#
my $reduceFile =$ARGV[0];
my $containerNumber=0;
my %containers;
open(FILE,$reduceFile);
while($line=<FILE>){
  chomp($line);
  my ($containee, $container, $ori, $offset)=split(/\s+/,$line);
  if(not(defined($containers{$container}))){
    $containers{$container}=$containerNumber;
    $containerNumber++;
  }
  $reduced{$containee}=$containers{$container}." $ori $offset";
}
close($reduceFile);
#we are only interested in containers/containees, the others are irreducible and thus cannot be contained
while($line=<STDIN>){
  chomp($line);
  my ($rname,$srname,$offset,$ori,$status)=split(/\s+/,$line);
  if(defined($containers{$srname})){#if container we can save space by renaming it to a number
    print "$rname ",$containers{$srname}," $offset $ori\n";
  }elsif(defined($reduced{$srname})){
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
