#!/usr/bin/env perl
#this code recovers the original scaffolds after gap closing -- all gaps are 100N's
my $ctgName="";
my $scf="";
my $chunk="";
my $flip=0;
while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
    $ctgName=substr($line,1);
    my @f=split(/\./,$ctgName);
    $scf=$f[0];
    #$scf=join(".",@f[0..($#f-1)]) if($#f>1);
    if($#f>1){#joined chunk, need to find orientation
      my @fb=split(/:/,$f[1]);
      my @fe=split(/:/,$f[-1]);
      $flip=1 if($fb[0]>$fe[0]);
    }
    my @ff=split(/:/,$f[-1]);
    $chunk=$ff[0];   
    $scfChunks{$scf}.="$chunk ";
  }else{
    if($flip){
      $line=reverse($line);
      $line=~tr/acgtACGT/tgcaTGCA/;
      $flip=0;
    }
    $ctgSeq{$scf.".".$chunk}=$line;
  }
}

#now we output the scaffolds
foreach $scf(keys %scfChunks){
my @f=split(/\s+/,$scfChunks{$scf});
if($#f==0){#only one chunk
  print ">$scf\n";
  print $ctgSeq{$scf.".".$f[0]},"\n";
}else{
  my @sorted_chunks=sort {$a <=> $b} @f;
  print ">$scf\n";
  print $ctgSeq{$scf.".".$sorted_chunks[0]};
    for($i=1;$i<=$#sorted_chunks;$i++){
      print "N"x100;
      print $ctgSeq{$scf.".".$sorted_chunks[$i]};
    }
    print "\n";
  }
}

