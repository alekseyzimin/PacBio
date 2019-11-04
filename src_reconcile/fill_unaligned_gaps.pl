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
my $mingap=1000; #minimum gap to be filled
my $gapnum=0;
my $gapbeg=0;
while($line=<STDIN>){
chomp($line);
@f=split(/\s+/,$line);
if($f[3]<$f[4]){
  $gapbeg=$f[0]-$f[3]+1;
}else{
  $gapbeg=$f[0]-($f[12]-$f[3]);
}
if($f[-2] eq $prevref && $gapbeg-$prevend>$mingap){#we found a fillable gap
  my $filllen=$gapbeg-$prevend-1;
  die("reference $f[-2] not found") if(not(defined($rseq{$f[-2]})));
  print STDERR ">fill$gapnum\n",substr($rseq{$f[-2]},$prevend,$filllen),"\n";
  print $prevend+1," ",$gapbeg-1," | 1 $filllen | $filllen $filllen | 100.0 | $f[11] $filllen | .1 100.0 | $f[-2] fill$gapnum\n";
  $gapnum++;
}
$prevref=$f[-2];
if($f[3]<$f[4]){
  $prevend=$f[1]+($f[12]-$f[4]);
}else{
  $prevend=$f[1]+$f[4]-1;
}
print $line,"\n";
}

