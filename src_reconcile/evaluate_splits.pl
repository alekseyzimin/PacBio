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
  @breaks=();
  @lines=();
  while($line=<STDIN>){
    chomp($line);
    my $ctg="";
    if($line=~/^break/ || $line=~/^alnbreak/){
      push(@breaks,$line);
    }elsif($line=~/^\-\-/){
      my $mincov=1000;
      my $maxcov=0;
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
        if($f[-1]>$maxcov){
          $maxcov=$f[-1];
        }
      }
      if(not($mincovline eq "")){
        @f=split(/\s+/,$mincovline);
        if($f[2]<5000 ||$f[2]>$sizes{$f[1]}-5000){
          print "$f[0] end_too_close_$f[1] $f[2] $f[3] ",join(" ",@breaks),"\n";
        }else{
          print "$mincovline ",join(" ",@breaks),"\n";
        }
        }
        foreach $b(@breaks){
          my @ff=split(/\s+/,$b);
          if($f[2]<5000 ||$f[2]>$sizes{$f[1]}-5000){
            print "repeat end_too_close_$ff[1] $ff[2] $maxcov ",join(" ",@breaks),"\n";
          }else{
            print "repeat $ff[1] $ff[2] $maxcov ",join(" ",@breaks),"\n";
          }
        }
        last;
    }else{
      push(@lines,$line);
    }
  }
}
