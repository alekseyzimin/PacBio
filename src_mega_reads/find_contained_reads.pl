#!/usr/bin/env perl
#
#this code, given the output of the super-mega-reads, eliminate contained mega-reads
my $readPlacementFile=$ARGV[0];
my $megaReadsFile=$ARGV[1];
my @mr_sizes;
my @mr_names;
my %groups;

open(FILE,$megaReadsFile);
while($line=<FILE>){
if(substr($line,0,1) eq ">"){
chomp($line);
push(@mr_names,substr($line,1));
}else{
push(@mr_sizes,length($line));
}
}

open(FILE,$readPlacementFile);
while($line=<FILE>){
chomp($line);
@f=split(/\s+/,$line);
my $mrn=int(substr($f[0],2)/2);
if($f[3] eq "F"){
$groups{$f[1]}.="$mrn ".($f[2]+1)." ".($mr_sizes[$mrn]+$f[2]-1)." ";
}else{
$groups{$f[1]}.="$mrn ".($f[2]-$mr_sizes[$mrn]+2)." $f[2] ";
}
}


die("error reading mega-reads file") if(not($#mr_sizes==$#mr_names));

foreach $g(keys %groups){
@f=split(/\s+/,$groups{$g});
$max_len=0;
$max_name="";
for($i=0;$i<$#f;$i+=3){
#do all vs all within each group
for($j=0;$j<$#f;$j+=3){
next if($i==$j);
next if($contained{$f[$j]});
next if($f[$j+2]-$f[$j+1]>$f[$i+2]-$f[$i+1]);#next if ith read is shorter than the jth read
if($f[$j+1]>=$f[$i+1] && $f[$j+2]<=$f[$i+2]){
print "$mr_names[$f[$j]]\n";
$contained{$f[$j]}=1;
}
}
}
}
