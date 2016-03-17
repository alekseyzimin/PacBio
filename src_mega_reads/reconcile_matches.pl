#!/usr/bin/env perl


my $tol_factor=1.2;#tolerance factor around gap size
my $tol_min=100;#minimum tolerance for small gap
#gap is an actual gap; all coordinates in the output are 1-based
my $scf="";
my $line=<>,@l=();
my @f=split(/\s+/,$line);
my $scf=$f[-2];
push(@l,$line);

while($line=<>){
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
  my $gap_before=10000;
  my $gap_after=10000;
  if($#lines==0){
     @l2=split(/\s+/,$lines[0]);
    output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
  }else{
    my @l1,@l2,@l3;
    @l2=split(/\s+/,$lines[0]);
    @l3=split(/\s+/,$lines[1]);
    $gap_after=compute_gap($l2[1],$l3[0]);
    output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
    for($i=1;$i<$#lines;$i++){
      @l1=split(/\s+/,$lines[$i-1]);
      @l2=split(/\s+/,$lines[$i]);
      @l3=split(/\s+/,$lines[$i+1]);

      $gap_before=compute_gap($l1[1],$l2[0]);
      $gap_after=compute_gap($l2[1],$l3[0]);

      output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
    }

    @l1=split(/\s+/,$lines[$#lines-1]);
    @l2=split(/\s+/,$lines[$#lines]);
    $gap_before=compute_gap($l1[1],$l2[0]);
    $gap_after=10000;
    output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],0);
  }
}

sub output_coords{
  my $start,$end,$dir;
  my ($gap_b,$gap_a,$s,$e,$len,$scf,$ctg)=@_;
  if($s<$e){#forward match
    $dir="f";
    if($s-1<=$gap_b){
      $start=1;
      $gap_b-=($s-1);
      $gap_b=10 if($gap_b<1);
    }else{
      $start=$s;
    }
    if(($len-$e)<=$gap_a){
      $end=$len;
      $gap_a-=($len-$e);
      $gap_a=10 if($gap_a<1);
    }else{
      $end=$e;
    }
  }else{
    $dir="r";
    if($e-1<=$gap_a){
      $start=1;
      $gap_a-=($e-1);
      $gap_a=10 if($gap_a<1);
    }else{
      $start=$e;
    }
    if(($len-$s)<=$gap_b){
      $end=$len;
      $gap_b-=($len-$s);
      $gap_b=10 if($gap_b<1);
    }else{
      $end=$s;
    }
  }
  print "$scf $ctg $start $end $dir ",int($gap_b)," ",int($gap_a)," $len\n";
}

sub compute_gap{
  my($gbeg,$gend)=@_;#gbeg<gend normally
    if($gend-$gbeg<1){
      return(1000);
    }elsif($gend-$gbeg==1){
      return(0);
    }else{
      my $ttt=($gend-$gbeg)*$tol_factor;
      return($ttt<$tol_min ? $tol_min : $ttt);
    }
}

