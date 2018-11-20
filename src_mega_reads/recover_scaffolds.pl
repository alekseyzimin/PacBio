#!/usr/bin/env perl
#this code recovers the original scaffolds after gap closing -- all gaps are 100N's
while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
    $ctgName=substr($line,1);
    my ($scf,$chunk)=split(/\./,$ctgName);
    $scfChunks{$scf}.="$chunk ";
  }else{
    $ctgSeq{$ctgName}=$line;
  }
}

#now we output the scaffolds
foreach $scf(keys %scfChunks){
my @f=split(/\s+/,$scfChunks{$scf});
if($#f==0){#only one chunk
  print ">$scf\n";
  print $ctgSeq{$scf.".".$f[0]},"\n";
}else{
  my @sorted_chunks=sort by_chunk(@f);
  print ">$scf\n";
  print $ctgSeq{$scf.".".$sorted_chunks[0]};
    for($i=1;$i<=$#sorted_chunks;$i++){
      print "N"x100;
      print $ctgSeq{$scf.".".$sorted_chunks[$i]};
    }
    print "\n";
  }
}

sub by_chunk{
  my @f1=split(/:/,$a);
  my @f2=split(/:/,$b);
  return($f1[0]<=>$f2[0]);
}
