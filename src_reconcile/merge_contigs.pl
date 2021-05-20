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
    my($ctg1,$dir1,$ctg2,$dir2,$gap)=split(/\s+/,$line);
    next if($gap>$max_gap);
    if($dir1 eq "F"){
        $edge_fwd{$ctg1}= defined($edge_fwd{$ctg1}) ? -1 : "$ctg2 $dir2 $gap";	
	if($dir2 eq "F"){
	    $edge_rev{$ctg2}=defined($edge_rev{$ctg2}) ? -1 : "$ctg1 F $gap";
	}else{
	    $edge_fwd{$ctg2}=defined($edge_fwd{$ctg2}) ? -1 : "$ctg1 R $gap";
	}
    }else{
	my $tdir=($dir2 eq "F") ? "R" : "F";
	$edge_rev{$ctg1}=defined($edge_rev{$ctg1}) ? -1 : "$ctg2 $tdir $gap";
	if($dir2 eq "F"){
	    $edge_rev{$ctg2}=defined($edge_rev{$ctg2}) ? -1 : "$ctg1 R $gap";
	}else{
	    $edge_fwd{$ctg2}=defined($edge_fwd{$ctg2}) ? -1 : "$ctg1 F $gap";
	}
    }
}
#now we delete all hash entries with -1
my @temp=keys %edge_fwd;
foreach my $e(@temp){
  delete $edge_fwd{$e} if($edge_fwd{$e}==-1);
#  print "fwd $e $edge_fwd{$e}\n" if(defined($edge_fwd{$e}));
}
my @temp=keys %edge_rev;
foreach my $e(@temp){
  delete $edge_rev{$e} if($edge_rev{$e}==-1);
#  print "rev $e $edge_rev{$e}\n" if defined($edge_rev{$e});;
}
#and delete all non-reciprocal edges
my @temp=keys %edge_fwd;
foreach my $e(@temp){
  my ($c,$d,$g)=split(/\s+/,$edge_fwd{$e});
  if($d eq "F"){
    delete  $edge_fwd{$e} if not(defined($edge_rev{$c}));
  }else{
    delete  $edge_fwd{$e} if not(defined($edge_fwd{$c}));
  }
}
my @temp=keys %edge_rev;
foreach my $e(@temp){
  my ($c,$d,$g)=split(/\s+/,$edge_rev{$e});
  if($d eq "F"){
    delete  $edge_rev{$e} if not(defined($edge_fwd{$c}));
  }else{
    delete  $edge_rev{$e} if not(defined($edge_rev{$c}));
  }
} 



#traverse the graph, go though edges, find terminal nodes and construct the graph
my $path="";
foreach my $e(keys %edge_fwd){
  next if(defined($edge_rev{$e})); #skip if internal node
  next if(defined($ctg_used{$e}));
  $ctg_used{$e}=1;
  $path="$e F ";
#print "DEBUG $path\n";
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
    $path.="$g $c $d ";
#print "DEBUG $path\n";
    $current_dir=$d;
die("fork detected in the forward loop $c |$path") if($ctg_used{$c});
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
die("fork detected in the reverse loop $c |$path") if($ctg_used{$c});
    $ctg_used{$c}=1;
  }while(defined($edge_rev{$c}) && defined($edge_fwd{$c}) && $last==0);
  $path=~s/^\s//;
  print $path,"\n";
}

