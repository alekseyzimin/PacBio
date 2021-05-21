#!/usr/bin/env perl
#this code takes contig graph in the format
#c1 dir c2 dir gap
#and then produces strings of merged contigs
#
my $max_gap=10000000;
my %edge_fwd;
my %edge_rev;
my %ctg_used;
my @links;
my $contigs=$ARGV[0];# we need contig sizes for bubble popping

open(FILE,$ARGV[0]);
my $len=-1;
while($line=<FILE>){
  chomp($line);
  if($line=~/^>/){
    my @f=split(/\s+/,$line);
    $len{$ctg}=$len if($len>-1);
    $ctg=substr($f[0],1);
    $len=0;
  }else{
    $len+=length($line);
  }
}
$len{$ctg}=$len if($len>-1);

#read in the graph
while($line=<STDIN>){
    chomp($line);
    my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$line);
    next if($gap>$max_gap);
    push(@links, $line);#save the links
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
#the idea is to first detect and collapse all linear paths
#so we delete all branches
#now we delete all hash entries with -1
my @temp=keys %edge_fwd;
foreach my $e(@temp){
  delete $edge_fwd{$e} if($edge_fwd{$e}==-1);
}
my @temp=keys %edge_rev;
foreach my $e(@temp){
  delete $edge_rev{$e} if($edge_rev{$e}==-1);
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

#find linear paths in the graph
my $path="";
my $pathindex=0;
my %path_beg=(), %path_end=(),@paths=();
walk_graph();

#now we create new links where the contigs are the paths
for($i=0;$i<=$#links;$i++){
  my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$links[$i]);
  my $origline="$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap";
  #print "reexamine $origline\n";
  my $tdir1=($dir1 eq "F") ? "R" : "F";
  my $tdir2=($dir2 eq "F") ? "R" : "F";

  if(defined($path_end{$ctg1.$dir1})){
    $ctg1="path".$path_end{$ctg1.$dir1};
    $dir1="F";
  }elsif(defined($path_beg{$ctg1.$tdir1})){
    $ctg1="path".$path_beg{$ctg1.$tdir1};
    $dir1="R";
  }
  if(defined($path_beg{$ctg2.$dir2})){
    $ctg2="path".$path_beg{$ctg2.$dir2};
    $dir2="F";
  }elsif(defined($path_end{$ctg2.$tdir2})){
    $ctg2="path".$path_end{$ctg2.$tdir2};
    $dir2="R";
  }
  $links[$i]="$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap" if(not("$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap" eq $origline));
  #print "NEW $ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap\n" if(not("$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap" eq $origline));
}

#now we go through the links again and look for nodes that form simple bubbles
#simple bubble node must have one forward and one reverse link
#forward and reverse links mus be to the same two in-degree 2 nodes

foreach $line(@links){
  my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$line);
  if($dir1 eq "F"){
    my $tdir=($dir2 eq "F") ? "R" : "F";
    next if defined $edge_fwd{$ctg1};#linear path edge, already used
    $edge_fwd_b{$ctg1}.="$ctg2 $dir2 $gap ";
    if($dir2 eq "F"){
      $edge_rev_b{$ctg2}.="$ctg1 F $gap ";
    }else{
      $edge_fwd_b{$ctg2}.="$ctg1 R $gap";
    }
  }else{
    my $tdir=($dir2 eq "F") ? "R" : "F";
    next if defined $edge_rev{$ctg1};#linear path edge, already used
    $edge_rev_b{$ctg1}.="$ctg2 $tdir $gap ";
    if($dir2 eq "F"){
      $edge_rev_b{$ctg2}="$ctg1 R $gap ";
    }else{
      $edge_fwd_b{$ctg2}="$ctg1 F $gap ";
    }
  }
}
#now edge_rev_b and edge_fwd_b contain all links with multiple connections
#we are looking for nodes that have exactly one forward and one reverse link
foreach $c(keys %edge_fwd_b){
  next if not defined $edge_rev_b{$c};
  my @fwd=split(/\s+/,$edge_fwd_b{$c});
  my @rev=split(/\s+/,$edge_rev_b{$c});
  $bubble{"$rev[0] $rev[1] $fwd[0] $fwd[1]"}.="$c " if($#fwd==2 && $#rev==2 && not($fwd[0] eq $rev[0]));
  $bubbleinfo{$c}="$edge_rev_b{$c} $edge_fwd_b{$c}";
}

foreach $k(keys %bubble){
  my @f=split(/\s+/,$bubble{$k});
  if($#f>0){
    my $bctg;
#pop the bubble
    if($len{$f[0]}>$len{$f[1]}){
      print STDERR "$f[1]\n";
      $bctg=$f[0];
    }else{
      print STDERR "$f[0]\n";
      $bctg=$f[1];
    }
    #print "bubble $k $bubble{$k}\n";
    my($ctg1,$dir1,$gap1,$ctg2,$dir2,$gap2)=split(/\s+/,$bubbleinfo{$bctg});
#now we record the four additional edges for each bubble
    #print "$ctg1,$dir1,$gap1,$ctg2,$dir2,$gap2\n";
    if($ctg1 =~ /^path/){
      my $pathindex=substr($ctg1,4);
      my @fff=split(/\s+/,$paths[$pathindex]);
      if($dir1 eq "F"){#end of the path
        $ctg1=$fff[-2];
        $dir1=$fff[-1];
      }else{
        $ctg1=$fff[0];
        $dir1=($fff[1] eq "F") ? "R" : "F";
      }
    }
    if($ctg2 =~ /^path/){
      my $pathindex=substr($ctg2,4);
      my @fff=split(/\s+/,$paths[$pathindex]);
      if($dir2 eq "F"){#begin of the path
        $ctg2=$fff[0];
        $dir2=$fff[1];
      }else{
        $ctg2=$fff[-2];
        $dir2=($fff[-1] eq "F") ? "R" : "F";
      }
    }
    #print "recover $bctg $ctg1,$dir1,$gap1,$ctg2,$dir2,$gap2\n";
    $edge_rev{$bctg}="$ctg1 $dir1 $gap1 ";
    $edge_fwd{$bctg}="$ctg2 $dir2 $gap2 ";
#reciprocals
    if($dir2 eq "F"){
      $edge_rev{$ctg2}="$bctg F $gap2 ";
    }else{
      $edge_fwd{$ctg2}="$bctg R $gap2 ";
    }
    if($dir1 eq "F"){
      $edge_fwd{$ctg1}="$bctg F $gap1 ";
    }else{
      $edge_rev{$ctg1}="$bctg R $gap1 ";
    }
  }
}

#re-walk the graph
walk_graph();

#output the final paths
foreach $p(@paths){
  print "$p\n";
}

#this sub walks the graph assuming all paths are linear
sub walk_graph{
$path="";
$pathindex=0;
%path_beg=(), %path_end=(),@paths=(),%ctg_used=();
#look at forward starting edges
foreach my $e(keys %edge_fwd){
  next if(defined($edge_rev{$e})); #skip if internal node
  next if(defined($ctg_used{$e}));
  $ctg_used{$e}=1;
  $path="$e"." F ";
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
    $current_dir=$d;
die("fork detected in the forward loop $c |$path") if($ctg_used{$c});
    $ctg_used{$c}=1;
  }while(defined($edge_rev{$c}) && defined($edge_fwd{$c}) && $last==0);
  push(@paths,$path);
  my @f=split(/\s+/,$path);
  $path_beg{$f[0].$f[1]}=$pathindex;
  $path_end{$f[-2].$f[-1]}=$pathindex;
  $pathindex++;
}
#look at th remaining reverse starting edges
foreach my $e(keys %edge_rev){
  next if(defined($edge_fwd{$e})); #skip if internal node
  next if(defined($ctg_used{$e}));
  $ctg_used{$e}=1;
  $path=" $e "."F";
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
  push(@paths,$path);
  my @f=split(/\s+/,$path);
  $path_beg{$f[0].$f[1]}=$pathindex;
  $path_end{$f[-2].$f[-1]}=$pathindex;
  $pathindex++;
}
}

