#!/usr/bin/env perl
#ths code fills unaligned gaps with reference sequence
#
my $refseq=$ARGV[0];

open(FILE,$refseq);
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

my $prevref;
my $prevend;
my $mingap=500;
my $gapnum=0;
while($line=<STDIN>){
chomp($line);
@f=split(/\s+/,$line);
if($f[-2] eq $prevref && $f[0]-$prevend>$mingap){#we found a fillable gap
my $filllen=$f[0]-$prevend-1;
die("reference $f[-2] not found") if(not(defined($rseq{$f[-2]})));
print STDERR ">fill$gapnum\n",substr($rseq{$f[-2]},$prevend+1,$filllen),"\n";
print $prevend+1," ",$f[0]-1," | 1 $filllen | $filllen $filllen | 100.0 | $f[11] $filllen | .1 100.0 | $f[-2] fill$gapnum\n";
$gapnum++;
}
$prevref=$f[-2];
$prevend=$f[1];
print $line,"\n";
}

