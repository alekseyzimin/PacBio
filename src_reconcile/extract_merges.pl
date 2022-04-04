#!/usr/bin/env perl
#
#this code extracts possible contig merges from a nucmer alignment: the reference sequences are merged with query sequence intput is a delta file 
#ASSUMES show-coords -q output, that is sorted by query coord!!!!!!

my $maxgap=500000;
my $mingap=-10000;

my $output_prefix="patches";
open(FILE,$ARGV[0]);#file with query contigs
while($line=<FILE>){
  chomp($line);
  if($line =~ /^>/){
    my @f=split(/\s+/,$line);
    $qn=substr($f[0],1);
  }else{
    $qseq{$qn}.=$line;
  }
}

my $min_match=500;
if(defined($ARGV[1])){
  $min_match=$ARGV[1];
}
my $max_overhang=1000;
if(defined($ARGV[2])){
  $max_overhang=$ARGV[2];
}

my $type="ont";
if(defined($ARGV[3])){
  $type=$ARGV[3];
} 

my $only_allowed=0;
my %allowed_merges=();
if(defined($ARGV[4])){
  $only_allowed=1;
  $maxgap=50*$max_overhang;
  $mingap=-1*$max_overhang;
  open(FILE,$ARGV[4]);
  while($line=<FILE>){
    chomp($line);
    my @f=split(/\s+/,$line);
    $allowed{"$f[0] $f[1]"}=1 if($#f==1);
  }
}

#first we read in all the matches into an array
my $prevline="";
my %oh1=(),%oh2=(),%gap=(),%gapseq=(),%idy=();

while($line=<STDIN>){
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  #this is the location where we filter the alignments by match length
  #we also allow to merge with a contig where over 95% of its length is aligned
  #filtering by overhang is provided earlier
  #we use matches that are smaller than the min_match to collect reads that span the gap, but may not pass the minimum match criteria, for consensus
  push(@lines,$line) if(($f[7] >= $min_match || $f[14]>95) && ($f[0]-1<=$max_overhang || $f[1]>=$f[11]-$max_overhang));
  $read_on_contig{$f[-2]}.="$f[-1] ";#for each contig get a list of reads that map to it
}

#now we go through the array and collect the possible merges
#
#
my $max_offset=1;
$max_offset=30 if($only_allowed);
for($i=0;$i<$#lines;$i++){
  @f1=split(/\s+/,$lines[$i]);
  for($j=$i+1;$j<=$i+$max_offset;$j++){
    next if($j>$#lines);
#print "DEBUG $i $j\n$lines[$i]\n$lines[$j]\n"; 
    @f2=split(/\s+/,$lines[$j]);
    next if($f1[-2] eq $f2[-2]);
    next if(not($f1[-1] eq $f2[-1]));
    my $gstart=1;
    my $dir1,$dir2,$oh1,$oh2;
#print "DEBUG considering $i $j\n$lines[$i]\n$lines[$j]\n\n";
    if($f1[3]<$f1[4]){
      $gstart=$f1[4];
      if($f2[3]<$f2[4]){
#forward forward merge ---->     ------>
        $gap=$f2[3]-$f1[4]+1;
        $oh1=$f1[11]-$f1[1];
        $oh2=$f2[0]-1;
        $dir1="F";
        $dir2="F";
      }else{
#forward reverse merge ---->     <------
        $gap=$f2[4]-$f1[4]+1;
        $oh1=$f1[11]-$f1[1];
        $oh2=$f2[11]-$f2[1];
        $dir1="F";
        $dir2="R";
      }
    }else{
      $gstart=$f1[3];
      if($f2[3]<$f2[4]){
#reverse forward merge <-----     ------>
        $gap=$f2[3]-$f1[3]+1;
        $oh1=$f1[0]-1;
        $oh2=$f2[0]-1;
        $dir1="R";
        $dir2="F";
      }else{
#reverse reverse merge <-----     <------
        $gap=$f2[4]-$f1[3]+1;
        $oh1=$f1[0]-1;
        $oh2=$f2[11]-$f2[1];
        $dir1="R";
        $dir2="R";
        } 
    }
    if($only_allowed){
      #check if the merge is allowed
      #allow for a single contig to be flipped
      next if((not(defined($allowed{"$f1[-2] $f2[-2]"})) && not(defined($allowed{"$f2[-2] $f1[-2]"}))) || (defined($allowed{"$f1[-2] $f2[-2]"}) && $dir1 eq "R" && $dir2 eq "R") || (defined($allowed{"$f2[-2] $f1[-2]"}) && $dir1 eq "F" && $dir2 eq "F"));
    }
    #print "DEBUG PASS $gap $oh1 $oh2\n";
    if($gap < $maxgap && $gap > $mingap && $oh1<=$max_overhang && $oh2<=$max_overhang){
        #$j=$i+$max_offset if($dir1 eq $dir2);#stop skipping search if found a candidate merge in the correct orientation for closing gaps
        $gstart=1 if($gstart<1);
        my $jstart=0;
        my $jend=length($qseq{$f1[-1]});
        if($type eq "asm"){
          my $fudge=5;
          $jstart = $gstart-1-$min_match*$fudge-$max_overhang > 0 ? $gstart-1-$min_match*$fudge-$max_overhang : 0;
          if($gap >= 0){
            $jend = $gstart-1+$gap+$min_match*$fudge+$max_overhang <= length($qseq{$f1[-1]}) ? $gstart-1+$gap+$min_match*$fudge+$max_overhang : length($qseq{$f1[-1]});
          }else{
            $jend = $gstart-1+$min_match*$fudge+$max_overhang <= length($qseq{$f1[-1]}) ? $gstart-1+$min_match*$fudge+$max_overhang : length($qseq{$f1[-1]});
          }
        }
        #print "DEBUG $jstart $jend\n";
        if($f1[-2] lt $f2[-2]){
          $joinline="$f1[-2]:$dir1:$f2[-2]:$dir2";
          #print "DEBUG joinline $joinline $oh1 $oh2 $gap\n";
          if(not(defined($oh1{$joinline})) || $oh1{$joinline}+$oh2{$joinline} > $oh1+$oh2 ){#here we use the best join for each pair of contigs, we minimize sum of overhangs
            $gseq{$joinline}= $gap>0 ? lc(substr($qseq{$f1[-1]},$gstart-1,$gap)) : "n";
            $jseq{$joinline}=lc(substr($qseq{$f1[-1]},$jstart,$jend-$jstart));
            $oh1{$joinline}=$oh1;
            $oh2{$joinline}=$oh2;
            $gap{$joinline}=$gap;
          }
          $paircount{"$f1[-2] $f2[-2]"}++;
        }else{
          $dir1= $dir1 eq "F" ? "R" : "F";
          $dir2= $dir2 eq "F" ? "R" : "F";
          $joinline="$f2[-2]:$dir2:$f1[-2]:$dir1";
          #print "DEBUG joinline $joinline $oh1 $oh2 $gap\n";
          if(not(defined($oh1{$joinline})) || $oh1{$joinline}+$oh2{$joinline} > $oh1+$oh2 ){#here we use the best join for each pair of contigs, we minimize sum of overhangs
            $gseq{$joinline}= $gap>0 ? reverse_complement(lc(substr($qseq{$f1[-1]},$gstart-1,$gap))) : "n";
            $jseq{$joinline}=lc(substr($qseq{$f1[-1]},$jstart,$jend-$jstart));
            $oh1{$joinline}=$oh2;
            $oh2{$joinline}=$oh1;
            $gap{$joinline}=$gap;
          }
          $paircount{"$f2[-2] $f1[-2]"}++;
        }
      $joincount{$joinline}++;
      $rnames{$joinline}.="$f1[-1] ";
    }
  }#$j loop
}#$i loop

#here we remove keys from rnames if "only allowed" and there are orientation conflicts
if($only_allowed){
  my %fwd=();
  my %rev=();
  my @to_delete=();
  foreach my $k (keys %rnames){
    my @f=split(/:/,$k);
    if($f[1] eq $f[3]){
      $fwd{"$f[0] $f[2]"}=1;
    }else{
      $rev{"$f[0] $f[2]"}=1;
    }
  }
  foreach $k (keys %rnames){
    my @f=split(/:/,$k);
    #print "DEBUG checking $k\n";
    if(defined($fwd{"$f[0] $f[2]"}) && defined($rev{"$f[0] $f[2]"})){
      push(@to_delete,$k) if(not($f[1] eq $f[3]));
    }
  }
  foreach $k (@to_delete){
    #print "DEBUG deleted $k\n";
    delete($rnames{$k});
  }
}

my %jnames=();
#if asm then no consensus needed
if($type eq "asm"){
  if(-e "do_consensus.sh"){
    open(RAW,">patches.raw.fa");
    foreach my $jl(keys %jseq){#$jl is joinline, REF is the best joining sequence padded
      print RAW ">$jl\n$jseq{$jl}\n";
    }
  }
}else{
#now we have all information in the hashes, let's output the bundles
#the longest read is the seq, the rest used to polish
#we output the reads and run the consensus for each patch separately
  foreach my $k (keys %rnames){# $k is the joinline
    my @names=split(/\s+/,$rnames{$k});
    my $max_len=0;
    my $max_i=0;
    for(my $i=0;$i<=$#names;$i++){
      if(length($qseq{$names[$i]})>$max_len){
        $max_i=$i;
        $max_len=length($qseq{$names[$i]});
      }
    }
#here we found the longest joining sequence -- it will be our reference for the consensus
#note that it may be present more than once
#here we collect all non-redundant read names into %jnames hash
#and associate them with the longest read
#these can be done over multiple gaps, that is the same longest read can accumulate other reads from the other gaps
#%jnames will be a single entry of 1 if there is only one read in @names
    my %output=();
    $jnames{$names[$max_i]}=1 if(not(defined($jnames{$names[$max_i]})));
    $output{$names[$max_i]}=1;
    foreach my $n(@names){
      if(not(defined($output{$n}))){
        $jnames{$names[$max_i]}.=" $n";
        $output{$n}=1;
      }
    }
#here we also add all other reads that map to the two contigs
    my @f=split(/:/,$k);
    my @ff1=split(/\s+/,$read_on_contig{$f[0]});
    my @ff2=split(/\s+/,$read_on_contig{$f[2]});
    my %temp=();
    foreach my $n(@ff1){
      $temp{$n}=1;
    }
    foreach my $n(@ff2){
      if(not(defined($output{$n})) && defined($temp{$n})){
        $jnames{$names[$max_i]}.=" $n";
        $output{$n}=1;
      }
    }
  }
#output the seqs and run consensus
#do_consensus.sh must exist!
  if(-e "do_consensus.sh"){
    open(RAW,">patches.raw.fa");
    my $pindex=0;
    foreach my $name(keys %jnames){
      my @names=split(/\s+/,$jnames{$name});
      if($#names==0){#no polishing seq -- put into raw
        print RAW ">$name\n$qseq{$name}\n";
      }else{
        open(REF,">patches.ref.$pindex.fa");
        open(READS,">patches.reads.$pindex.fa");
        print REF ">$name\n$qseq{$name}\n";
        print READS ">_$name\n$qseq{$name}\n";
#need uniq to avoid outputting duplicates
        my %output=();
        for(my $i=1;$i<=$#names;$i++){
          if(not(defined($output{$names[$i]}))){
            print READS ">$names[$i]\n$qseq{$names[$i]}\n";
            print READS ">_$names[$i]\n$qseq{$names[$i]}\n" if($#names<5);
            $output{$names[$i]}=1;
          }
        }
        close(REF);
        close(READS);
        $pindex++;
#run consensus
        if($pindex>=10){
          system("./do_consensus.sh");
          $pindex=0;
        }
      }
    }
    if($pindex>0){
      system("./do_consensus.sh");
    }
  }
}

#output the links
foreach $k (keys %rnames){ #$k is the joinline
  my @f=split(/:/,$k);
  if($only_allowed){
    print "$f[0] $oh1{$k} $f[1] $f[2] $oh2{$k} $f[3] $gap{$k} $gseq{$k}\n";
  }elsif($paircount{"$f[0] $f[2]"} == $joincount{$k} || $joincount{$k}>1){
    print "$f[0] $oh1{$k} $f[1] $f[2] $oh2{$k} $f[3] $gap{$k} $gseq{$k}\n";
  }
}

sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}
