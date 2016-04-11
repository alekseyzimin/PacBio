#!/usr/bin/env perl


my $tol_factor=100000;#tolerance factor around gap size
my $tol_min=100;#minimum tolerance for small gap
#gap is an actual gap; all coordinates in the output are 1-based

my $gap_coordinates=$ARGV[0];
open(FILE,$gap_coordinates);

while(my $line=<FILE>){
chomp($line);
my @f=split(/\s+/,$line);
push(@{$gaps{$f[0]}},[$f[1],$f[2]]);
}

my $split_contigs=$ARGV[1]; #these are contigs that become duplicted if too much tolerance is allowed -- they are likely misassembled
open(FILE,$split_contigs);

while(my $line=<FILE>){
chomp($line);
@f=split(/\s+/,$line);
$restrict_reverse_ext{"$f[0] $f[1] $f[4]"}=$f[2];
$restrict_forward_ext{"$f[0] $f[1] $f[4]"}=$f[3];
}


my $scf="";
my $line=<STDIN>,@l=();
my @f=split(/\s+/,$line);
my $scf=$f[-2];
push(@l,$line);

while($line=<STDIN>){
  $line=~s/^\s+//;
  @f=split(/\s+/,$line);
  if(not($f[-2] eq $scf)){
    process_lines(@l);
    @l=();
    $scf=$f[-2];
  }
  push(@l,$line);
}
process_lines(@l);

sub process_lines{
  my @lines=@_;
  my $gap_before=500000;
  my $gap_after=10000;
  if($#lines==0){
     @l2=split(/\s+/,$lines[0]);
    output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
  }else{
    my @l1,@l2,@l3;
    @l2=split(/\s+/,$lines[0]);
    @l3=split(/\s+/,$lines[1]);
    $gap_after=compute_gap($l2[1],$l3[0],$l2[-2]);
    output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
    for($i=1;$i<$#lines;$i++){
      @l1=split(/\s+/,$lines[$i-1]);
      @l2=split(/\s+/,$lines[$i]);
      @l3=split(/\s+/,$lines[$i+1]);

      $gap_before=compute_gap($l1[1],$l2[0],$l2[-2]);
      $gap_after=compute_gap($l2[1],$l3[0],$l2[-2]);

      output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
    }

    @l1=split(/\s+/,$lines[$#lines-1]);
    @l2=split(/\s+/,$lines[$#lines]);
    $gap_before=compute_gap($l1[1],$l2[0],$l2[-2]);
    $gap_after=500000;
    output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],0);
  }
}

sub output_coords{
  my $start,$end,$dir;
  my ($gap_b,$gap_a,$s,$e,$len,$scf,$ctg)=@_;
  my $sg_a=0,$sg_b=0;
   
  if($gap_b<0){
    $gap_b=-$gap_b;
    $sg_b=1;
  }
  if($gap_a<0){
    $gap_a=-$gap_a;
    $sg_a=1;
  }

  if($s<$e){#forward match
    $dir="f";
    $min_coord=1;
    $min_coord=int($s-$restrict_reverse_ext{"$s $e $ctg"}) if(defined($restrict_reverse_ext{"$s $e $ctg"}));
    $min_coord=1 if($min_coord<1);
    if($s-$min_coord<=$gap_b){
      $start=$min_coord;
      $gap_b-=($s-$min_coord);
    }else{
      $start=$s;
    }
    $max_coord=$len;
    $max_coord=int($e+$restrict_forward_ext{"$s $e $ctg"}) if(defined($restrict_forward_ext{"$s $e $ctg"}));
    $max_coord=$len if($max_coord>$len);
    if(($max_coord-$e)<=$gap_a){
      $end=$max_coord;
      $gap_a-=($max_coord-$e);
    }else{
      $end=$e;
    }
  }else{
    $dir="r";
    $min_coord=1;
    $min_coord=int($e-$restrict_reverse_ext{"$e $s $ctg"}) if(defined($restrict_reverse_ext{"$e $s $ctg"}));
    $min_coord=1 if($min_coord<1);
    if($e-$min_coord<=$gap_a){
      $start=$min_coord;
      $gap_a-=($e-$min_coord);
    }else{
      $start=$e;
    }
    $max_coord=$len;
    $max_coord=int($s+$restrict_forward_ext{"$e $s $ctg"}) if(defined($restrict_forward_ext{"$e $s $ctg"}));
    $max_coord=$len if($max_coord>$len);

    if(($max_coord-$s)<=$gap_b){
      $end=$max_coord;
      $gap_b-=($max_coord-$s);
    }else{
      $end=$s;
    }
  }

  $gap_a=$gap_a/$tol_factor if($sg_a);
  $gap_b=$gap_b/$tol_factor if($sg_b);

  print "$scf $ctg $start $end $dir ",int($gap_b)," ",int($gap_a)," $len\n";
}

sub compute_gap{
  my($gbeg,$gend,$name)=@_;#gbeg<gend normally
  #check if this corresponds to a sequence gap
  my $seq_gap=0;
  foreach my $g(@{$gaps{$name}}){
    if($$g[0]>=$gbeg && $$g[1]<=$gend){
      $seq_gap=1;
      last;
    }elsif($g[0]>$gend){
      last;
    }
  }
  if($seq_gap){
    #print "found sequence gap $gbeg,$gend,$name\n";
    my $ttt=($gend-$gbeg)*$tol_factor;
    $ttt=$tol_min*$tol_factor if($ttt<$tol_min);
    return(-$ttt); #negative gap is a flag for sequence gap
  }elsif($gend-$gbeg==1){
    return(0);
  }elsif($gend-$gbeg<=0){
    return(-100*$tol_factor);
  }else{
    my $ttt=($gend-$gbeg);
    return($ttt<$tol_min ? -$tol_min*$tol_factor : -$ttt*$tol_factor);
  }
}

