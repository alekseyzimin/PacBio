#!/usr/bin/env perl
#here we take paths on stdin and the repeats file on command line and insert repeat contigs into the paths
my $repeat_info=$ARGV[0];
my %rep_before=(),%rep_after=();
open(FILE,$repeat_info);
while($line=<FILE>){
  chomp($line);
  my @f=split(/\s+/,$line);
  $rep_ctg=$f[0];
  if($f[1] eq "fwd"){
    for($i=2;$i<$#f;$i+=3){
      $rep_before{$f[$i].$f[$i+1]}="$f[0] $f[$i+2]";
    }
  }else{
    for($i=2;$i<$#f;$i+=3){
      $rep_after{$f[$i].$f[$i+1]}="$f[0] $f[$i+2]";
    }
  }
}

#read the paths and insert repeats
while($line=<STDIN>){
  my $newpath="";
  chomp($line);
  @f=split(/\s+/,$line);
  $newpath.="$f[0] $f[1] ";
  for($i=3;$i<=$#f;$i+=3){
    my $tdir1=$f[$i-2] eq "F" ? "R" : "F";
    my $tdir2=$f[$i+1] eq "F" ? "R" : "F";
    #print "TEST $rep_after{$f[$i-3].$f[$i-2]} | $rep_before{$f[$i].$f[$i+1]} | $rep_before{$f[$i-3].$tdir1} | $rep_after{$f[$i].$tdir2}\n";
    if(defined($rep_after{$f[$i-3].$f[$i-2]}) && defined($rep_before{$f[$i].$f[$i+1]})){
      my ($ctga,$ga)=split(/\s+/,$rep_after{$f[$i-3].$f[$i-2]});
      my ($ctgb,$gb)=split(/\s+/,$rep_before{$f[$i].$f[$i+1]});
      if($ctga eq $ctgb){
        #inserting repeat
        $newpath.="$ga $ctga F $gb ";
        #print "DEBUG inserted $ga $ctga F $gb\n";
      }else{
        #no repeat
        $newpath.="$f[$i-1] ";
      }
    }elsif(defined($rep_before{$f[$i-3].$tdir1}) && defined($rep_after{$f[$i].$tdir2})){
      my ($ctga,$ga)=split(/\s+/,$rep_after{$f[$i].$tdir2});
      my ($ctgb,$gb)=split(/\s+/,$rep_before{$f[$i-3].$tdir1});
      if($ctga eq $ctgb){
        #inserting repeat reverse
        $newpath.="$gb $ctga R $ga ";
        #print "DEBUG inserted $ga $ctga R $gb\n";
      }else{
        $newpath.="$f[$i-1] ";
      }
    }else{
      $newpath.="$f[$i-1] ";
    }
    $newpath.="$f[$i] $f[$i+1] ";
  }
  print $newpath,"\n";
}

      
    
