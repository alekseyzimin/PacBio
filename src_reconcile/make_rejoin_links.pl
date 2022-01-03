#!/usr/bin/env perl
#This code makes rejoin links file from broken/merged scaffolds file.  input is the output of ufasta sizes -H on the scaffolds file
@lines=();
while($line=<STDIN>){
next if($line=~ /:/);#next if merged
chomp($line);
my @f=split(/\s+/,$line);
my @namearr=split(/\./,$f[0]);
next if($#namearr == 0);
next if(not($namearr[1] =~ /\d+/));
push(@lines,"$namearr[0] $namearr[1] $f[1]");
}

@lines = sort byscaffold_then_position @lines;

sub byscaffold_then_position {
my @fa = split(/\s+/,$a);
my @fb = split(/\s+/,$b);
return($fa[0] cmp $fb[0] || $fa[1] <=> $fb[1]);
}

my $prevname="";
my $prevcoord=-1;
my $prevoffset=-1;
foreach $l(@lines){
#print "$l\n";
my @f=split(/\s+/,$l);
  if($f[0] eq $prevname && $f[1] == $prevoffset){
    print "$f[0].$prevcoord 0 F $f[0].$f[1] 0 F 0 n\n";
  }
  $prevname=$f[0];
  $prevcoord=$f[1];
  $prevoffset=$f[1]+$f[2];
}


