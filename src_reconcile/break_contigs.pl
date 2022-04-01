#!/usr/bin/env perl
#
#
my $breaks_file=$ARGV[0];
open(FILE, $breaks_file);
while($line=<FILE>){
  chomp($line);
  @f=split(/\s+/,$line);
  my $break_coord=$f[2];
  push(@{$break_coords{$f[1]}},$break_coord);
}

while($line=<STDIN>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      if(defined($break_coords{$ctg})){
        my $offset=0;
        for(my $i=0;$i<@{$break_coords{$ctg}};$i++){
          #should not break within 5kb near a gap
          my $region=substr($seq,${$break_coords{$ctg}}[$i]-5000,10000);
          if(not($region =~ /N/ || $region =~ /n/)){ 
            print ">$ctg.$offset\n";
            print substr($seq,$offset,${$break_coords{$ctg}}[$i]-$offset),"\n";
            $offset=${$break_coords{$ctg}}[$i];
          }
        }
        if($offset>0){
          print ">$ctg.$offset\n";
          print substr($seq,$offset),"\n";
        }else{
          print ">$ctg\n$seq\n";
        }
        }else{
          print ">$ctg\n$seq\n";
        }
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}

if(not($seq eq "")){
  if(defined($break_coords{$ctg})){
    my $offset=0;
    for($i=0;$i<@{$break_coords{$ctg}};$i++){
      my $region=substr($seq,${$break_coords{$ctg}}[$i]-5000,10000);
      if(not($region =~ /N/ || $region =~ /n/)){
        print ">$ctg.$offset\n";
        print substr($seq,$offset,${$break_coords{$ctg}}[$i]-$offset),"\n";
        $offset=${$break_coords{$ctg}}[$i];
      }
    }
    if($offset>0){
      print ">$ctg.$offset\n";
      print substr($seq,$offset),"\n";
    }else{
      print ">$ctg\n$seq\n";
    }
  }else{
    print ">$ctg\n$seq\n";
  }
}

