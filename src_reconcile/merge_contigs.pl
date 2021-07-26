#!/usr/bin/env perl
#this code takes contig graph in the format
#c1 dir c2 dir gap
#and then produces strings of merged contigs
#
my $num_bubbles=0;
my %edge_fwd=();
my %edge_rev=();
my %ctg_used=();
my @links=(),@pathlinks=();
my $contigs=$ARGV[0];# we need contig sizes for bubble popping
my @edge_fwd_b=(),@edge_rev_b=(),%bubble=(),%bubbleinfo=();
my %path_beg=(), %path_end=(),@paths=(),%ctg_used=();
my $max_tip=10000;

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

#read in the links
while($line=<STDIN>){
    chomp($line);
    my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$line);
    push(@links, $line);#save the links
}

#find linear paths in the graph
walk_graph();
#pop bubbles
pop_bubbles();
#re-walk the graph
while($num_bubbles < keys %bubbleinfo){
  walk_graph();
  pop_bubbles();
  $num_bubbles=keys %bubbleinfo;
}

#output the final paths
foreach $p(@paths){
  print "$p\n";
}
foreach $c(keys %bubbleinfo){
  print STDERR "$c\n";
}

#this sub prines tips -- short contigs that do not continue but create a fork

#this sub pops bubbles
sub pop_bubbles{
#initialize
@edge_fwd_b=(),@edge_rev_b=(),%bubble=();
#now we go through the links again and look for nodes that form simple bubbles
#simple bubble node must have one forward and one reverse link
#forward and reverse links mus be to the same two in-degree 2 nodes
  foreach $line(@pathlinks){
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
  foreach my $c(keys %edge_fwd_b){
    next if not defined $edge_rev_b{$c};
    my @fwd=split(/\s+/,$edge_fwd_b{$c});
    my @rev=split(/\s+/,$edge_rev_b{$c});
    $bubble{"$rev[0] $rev[1] $fwd[0] $fwd[1]"}.="$c " if(not($fwd[0] eq $rev[0]));
  }

#now we pop the bubble by deleting the links to one of the bubble contigs (or paths) from the @links
  foreach my $k(keys %bubble){
    my @f=split(/\s+/,$bubble{$k});
    if($#f>0){
      my $bctg = $len{$f[0]} > $len{$f[1]} ? $f[1] : $f[0];
#pop the bubble
#print "bubble $k $bubble{$k}\n";
      $bubbleinfo{$bctg}=$k;#this contains the links to be deleted
    }
  }

  for(my $i=0;$i<=$#links;$i++){
    my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$links[$i]);
#first we check of the contigs are beginnings or ends of the paths
    my $tdir1=($dir1 eq "F") ? "R" : "F";
    my $tdir2=($dir2 eq "F") ? "R" : "F";
    $ctg1="path".$path_end{$ctg1.$dir1} if(defined($path_end{$ctg1.$dir1}));
    if(defined($path_beg{$ctg1.$tdir1})){
      $ctg1="path".$path_beg{$ctg1.$tdir1};
      $dir1=$tdir1;
    }
    $ctg2="path".$path_beg{$ctg2.$dir2} if(defined($path_beg{$ctg2.$dir2}));
    if(defined($path_end{$ctg2.$tdir2})){
      $ctg2="path".$path_end{$ctg2.$tdir2};
      $dir2=$tdir2;
    }
#here we delete links to the bubble from the @links
    $links[$i]="" if(defined($bubbleinfo{$ctg1})||defined($bubbleinfo{$ctg2}));  
  }
}


#this sub walks the graph loading the links from @links by looking for all linear paths in @links
#it creates modified links, @pathlinks where all linear paths are collapsed to a single node
#then pathlinks are used in popping bubbles
sub walk_graph{
  my $path="";
  my $pathindex=0;
  @pathlinks=();
  %path_beg=(), %path_end=(),@paths=(),%ctg_used=();
  %edge_fwd=(), %edge_rev=();
  foreach $line(@links){
    next if($line eq "");
    push (@pathlinks,$line);
    my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$line);
    if($dir1 eq "F"){
      $edge_fwd{$ctg1}.="$ctg2 $dir2 $gap ";
      if($dir2 eq "F"){
        $edge_rev{$ctg2}.="$ctg1 F $gap ";
      }else{
        $edge_fwd{$ctg2}.="$ctg1 R $gap ";
      }
    }else{
      my $tdir=($dir2 eq "F") ? "R" : "F";
      $edge_rev{$ctg1}.="$ctg2 $tdir $gap ";
      if($dir2 eq "F"){
        $edge_rev{$ctg2}.="$ctg1 R $gap ";
      }else{
        $edge_fwd{$ctg2}.="$ctg1 F $gap ";
      }
    }
  }
#the idea is to first detect and collapse all linear paths
#so we delete all branches if they do not branch to a tip, otherwise we delete the tip
  my @temp=keys %edge_fwd;
  foreach my $e(@temp){
    my @f=split(/\s+/,$edge_fwd{$e});
    if($#f>2){#more than one edge, check for tips
      my %tips=();
      for($i=0;$i<$#f;$i+=3){
        @ff=split(/\s+/,$edge_fwd{$f[$i]});
        @fr=split(/\s+/,$edge_rev{$f[$i]});
        if((($#ff==2 && $#fr<2) || ($#ff<2 && $#fr==2)) && $len{$f[$i]}<$max_tip){#found a tip
          delete $edge_fwd{$f[$i]} if($#ff==2 && $#fr<2);
          delete $edge_rev{$f[$i]} if($#ff<2 && $#fr==2);
          $tips{$i}=1;
        }
      }
      my $newedges="";
      my $num_newedges=0;
      for($i=0;$i<$#f;$i+=3){
        if(not(defined($tips{$i}))){
          $newedges.="$f[$i] ".$f[$i+1]." ".$f[$i+2]." ";
          $num_newedges++;
        }
      }
      if($num_newedges==1){
        $edge_fwd{$e}=$newedges;
      }else{
        delete $edge_fwd{$e};
      }
    }
  }
  my @temp=keys %edge_rev;
  foreach my $e(@temp){
    my @f=split(/\s+/,$edge_rev{$e});

    if($#f>2){#more than one edge, check for tips
      my %tips=();
      for($i=0;$i<$#f;$i+=3){
        @ff=split(/\s+/,$edge_fwd{$f[$i]});
        @fr=split(/\s+/,$edge_rev{$f[$i]});
        if((($#ff==2 && $#fr<2) || ($#ff<2 && $#fr==2)) && $len{$f[$i]}<$max_tip){#found a tip
          delete $edge_fwd{$f[$i]} if($#ff==2 && $#fr<2);
          delete $edge_rev{$f[$i]} if($#ff<2 && $#fr==2);
          $tips{$i}=1;
        }
      }
      my $newedges="";
      my $num_newedges=0;   
      for($i=0;$i<$#f;$i+=3){
        if(not(defined($tips{$i}))){
          $newedges.="$f[$i] ".$f[$i+1]." ".$f[$i+2]." ";
          $num_newedges++;
        }
      }
      if($num_newedges==1){
        $edge_rev{$e}=$newedges;
      }else{
        delete $edge_rev{$e};
      }
    }
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
#now we walk the linear paths
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
#now we create new links where the contigs are the paths for bubble detection
  for($i=0;$i<=$#pathlinks;$i++){
    my($ctg1,$oh1,$dir1,$ctg2,$oh2,$dir2,$gap)=split(/\s+/,$pathlinks[$i]);
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
    $pathlinks[$i]="$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap" if(not("$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap" eq $origline));
#print "NEW $ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap\n" if(not("$ctg1 $oh1 $dir1 $ctg2 $oh2 $dir2 $gap" eq $origline));
  }
}

