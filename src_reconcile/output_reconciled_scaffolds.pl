#!/usr/bin/perl
my $name="";
my $seq="";
my $seqfile=$ARGV[0];
open(FILE,$seqfile);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    $sequence{$name}=$seq if(not($name eq ""));
    $output{$name}=0;
    ($name)=split(/\s+/,substr($line,1));
    $seq="";
  }else{
    $seq.=$line;
  }
}
$sequence{$name}=$seq;

$name="";
$seq="";
my $gap=1000000;
while($line=<STDIN>){
#print $line;
chomp($line);
@f=split(/\s+/,$line);
$gap=$gap<$f[5] ? $gap : $f[5];
if(not($f[0] eq $name)){
  print ">$name.R\n$seq\n" if(not($name eq ""));
  $name=$f[0];
  $seq="";
}else{
  $seq.=("N"x$gap) if($gap>0);
  #print "$gap\n";
}

die("Sequence $f[1] not found") if(not(defined($sequence{$f[1]})));
my $offset=1;
$offset=$gap+1 if($gap<0);
$seq.=($f[4] eq "f") ? substr($sequence{$f[1]},$f[2]-$offset,$f[3]-$f[2]+1) : reverse_complement(substr($sequence{$f[1]},$f[2]-$offset,$f[3]-$f[2]+1));
#print "$f[1] ",$f[2]-$offset," ",$f[3]-$f[2]+1," $f[4]\n";
$output{$f[1]}=1;
$gap=$f[6];
}
print ">$name.R\n$seq\n";

foreach $k(keys %output){
print ">$k\n$sequence{$k}\n" if(not($output{$k}) && length($sequence{$k})>1000);
}
sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}

