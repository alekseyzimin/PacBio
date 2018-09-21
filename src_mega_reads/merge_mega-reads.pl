#!/usr/bin/env perl
#this code takes contig graph in the format
#c1 dir c2 dir gap
#and then produces strings of merged contigs
#
my $max_gap=10000000;
my %edge_fwd;
my %edge_rev;
my %ctg_used;

#read in the graph
while($line=<STDIN>){
    chomp($line);
    my($ctg1,$trim1,$dir1,$ctg2,$trim2,$dir2,$gap)=split(/\s+/,$line);
    next if($gap>$max_gap);
    if($dir1 eq "F" && $dir2 eq "F"){
	$edge_fwd{$ctg1}="$ctg2 F $gap";
	$edge_rev{$ctg2}="$ctg1 F $gap";
    }else{
	$edge_rev{$ctg1}="$ctg2 F $gap";
        $edge_fwd{$ctg2}="$ctg1 F $gap";
    }
}


#traverse the graph, go though edges, find terminal nodes and construct the graph
my $path="";
foreach my  $e(keys %edge_fwd){
  next if(defined($edge_rev{$e})); #skip if internal node
  next if(defined($ctg_used{$e}));
  $ctg_used{$e}=1;
  $path="$e F ";
  my $current_dir="F";
  my $c=$e;
  my $last=0;
  do{
    if($current_dir eq "F"){
      ($c,$d,$g)=split(/\s+/,$edge_fwd{$c});
    }else{
      ($c,$d,$g)=split(/\s+/,$edge_rev{$c});
      $d=~tr/FR/RF/;
    }
    $last=1 if($ctg_used{$c});# if found a fork
#print "DEBUG $path $g $c $d | $current_dir\n";
    $path.="$g $c $d ";
    $current_dir=$d;
#die("fork detected in the forward loop $c $path") if($ctg_used{$c});
    $ctg_used{$c}=1;
  }while(defined($edge_rev{$c}) && defined($edge_fwd{$c}) && $last==0);
  print $path,"\n";
}

#print "reverse loop\n";
foreach my $e(keys %edge_rev){
  next if(defined($edge_fwd{$e})); #skip if internal node
  next if(defined($ctg_used{$e}));
  $ctg_used{$e}=1;
  $path=" $e F";
  my $current_dir="F";
  my $c=$e;
  my $last=0;
  do{
    if($current_dir eq "F"){
      ($c,$d,$g)=split(/\s+/,$edge_rev{$c});
    }else{
      ($c,$d,$g)=split(/\s+/,$edge_fwd{$c});
      $d=~tr/FR/RF/;
    }
    $last=1 if($ctg_used{$c});
    $path=" $c $d $g".$path;
    $current_dir=$d;
#die("fork detected in the reverse loop $c $path") if($ctg_used{$c});
    $ctg_used{$c}=1;
  }while(defined($edge_rev{$c}) && defined($edge_fwd{$c}) && $last==0);
  $path=~s/^\s//;
  print $path,"\n";
}

