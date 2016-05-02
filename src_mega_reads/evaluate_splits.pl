#!/usr/bin/env perl
#
my $contig_sizes=$ARGV[0];
open(FILE,$contig_sizes);
while($line=<FILE>){
chomp($line);
my ($ctg,$size)=split(/\s+/,$line);
$sizes{$ctg}=$size;
}

my @breaks;
my @lines;
while($line=<STDIN>){
  @abreaks=();
  @breaks=();
  @lines=();
  while($line=<STDIN>){
    chomp($line);
    my $ctg="";
    if($line=~/^b/ || $line=~/^a/){
      push(@breaks,$line);
    }elsif($line=~/^\-\-/){
      my $mincov=1000;
      my $mincovline="";
      my @f=split(/\s+/,$breaks[0]);
      $ctg=$f[1];
      foreach my $l(@lines){
        my @f=split(/\s+/,$l);
        next if($f[-1]==0 || $f[-2]<1000 || not($f[1] eq $ctg));
        if($f[-1]<$mincov){
          $mincov=$f[-1];
          $mincovline=$l;
        }
      }
      if(not($mincovline eq "")){
        @f=split(/\s+/,$mincovline);
        if($f[2]<5000 ||$f[2]>$sizes{$f[1]}-5000){
          print "$f[0] end_too_close_$f[1] $f[2] $f[3]";
        }else{
          print $mincovline;
        }
        print " ",join(" ",@breaks),"\n";
        }
        last;
    }else{
      push(@lines,$line);
    }
  }
}
