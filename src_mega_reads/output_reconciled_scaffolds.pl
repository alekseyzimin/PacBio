#!/usr/bin/perl
my $name="";
my $seq="";
my $seqfile=$ARGV[0];
open(FILE,$seqfile);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    $sequence{$name}=$seq if(not($name eq ""));
    ($name)=split(/\s+/,substr($line,1));
    $seq="";
  }else{
    $seq.=$line;
  }
}

$name="";
$seq="";
my $gap=1000000;
while($line=<STDIN>){
chomp($line);
@f=split(/\s+/,$line);
$gap=$gap<$f[5] ? $gap : $f[5];
if(not($f[0] eq $name)){
  print ">$name\n$seq\n" if(not($name eq ""));
  $name=$f[0];
  $seq="";
}else{
  $seq.=("N"x$gap);
}
$gap=$f[6];
die("Sequence $f[1] not found") if(not(defined($sequence{$f[1]})));
$seq.=($f[4] eq "f") ? substr($sequence{$f[1]},$f[2],$f[3]-$f[2]) : reverse_complement(substr($sequence{$f[1]},$f[2],$f[3]-$f[2]));
}
print ">$name\n$seq\n";

sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}

