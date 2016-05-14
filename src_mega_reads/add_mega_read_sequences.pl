#!/usr/bin/env perl
#this code adds missing mega-read sequences, one argument, k-unitigs file
my $kufile=$ARGV[0];
my @f,$kuname, $line,$kusize=10000000,@kunitigs;
open(FILE,$kufile);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    @f=split(/\s+/,$line);
    $kuname=substr($f[0],1)*2;
  }else{
    $kusize=length($line) if($kusize>length($line));
    $kunitigs[$kuname]=$line;
    $line=~tr/ACGTNacgtn/TGCANtgcan/;
    $line=reverse($line);
    $kunitigs[$kuname+1]=$line;
  }
}
$kusize--;#need to know k-1
my @kuext;
$#kuext=$#kunitigs;
for(my $i=0;$i<=$#kunitigs;$i++){
$kuext[$i]=substr($kunitigs[$i],$kusize);
}
#now we read in the mega reads file and add in the sequences
while($line=<STDIN>){
  chomp($line);
  if($line=~/^\>/){ 
    print $line,"\n";
  }else{
    my @f=split(/\s+/,$line);
    my @mrname=split(/_/,$f[8]);
    my $ku,$seq="";
    $ku=substr($mrname[0],0,-1)*2;
    if($mrname[0] =~ /F$/){
      $seq=$kunitigs[$ku];
    }else{
      $seq=$kunitigs[$ku+1];
    }
    for(my $i=1;$i<=$#mrname; $i++){
    $ku=substr($mrname[$i],0,-1)*2;
    if($mrname[$i]=~/F$/){
        $seq.=$kuext[$ku];
      }else{
        $seq.=$kuext[$ku+1];
      }
    }
    print "$line $seq\n";
  }
}

