#!/usr/bin/env perl
#first we read in the coords file to assign contigs to the chromosomes
while($line=<STDIN>){
  chomp($line);
  $line=~s/^\s+//;
  my @f=split(/\s+/,$line);
  next if($f[7]<1000);
  my $impl_start;
  my $impl_end;
  my $dir;
  if($f[3]<$f[4]){#forward match
    $impl_start=$f[0]-$f[3];
    $impl_end=$f[1]+($f[12]-$f[4]);
    $dir="+";
  }else{
    $impl_start=$f[0]-($f[12]-$f[3]);
    $impl_end=$f[1]+$f[4];
    $dir="-";
  }
  next if($impl_start < -1000000 || $impl_end > $f[11]+1000000);
  my $midpoint=($impl_start+$impl_end)/2;
  $matches{$f[18]}.="$f[17] $f[7] $midpoint $dir ";
  $ref_chr_len{$f[17]}=$f[11] if(not(defined($ref_chr_len{$f[17]})));
  $ctg_len{$f[18]}=$f[12] if(not(defined($ctg_len{$f[18]})));
  #print "DEBUG $line\n$f[-2] $f[7] $impl_start $impl_end\n";
}

#now we assign contigs to the chromosomes
foreach $c(keys %matches){
  my @f=split(/\s+/,$matches{$c});
  my %temp=();
  for($i=0;$i<$#f;$i+=4){
    $temp{$f[$i]}+=$f[$i+1];
  }
  my $max_match=0;
  my $max_chrom="";
  foreach $ch(keys %temp){
    #print "DEBUG contig $c chromosome $ch bases $temp{$ch}\n";
    if($temp{$ch}>$max_match){
      $max_match=$temp{$ch};
      $max_chrom=$ch;
    }
  }
  $chrom{$c}=$max_chrom;
}

foreach $c(keys %matches){
  my @f=split(/\s+/,$matches{$c});
  #we compute the weighted mean position
  my $weight_sum=0;
  my $sum=0;
  my $count=0;
  my $fwd_sum=0;
  my $rev_sum=0;
  for($i=0;$i<$#f;$i+=4){
    if($chrom{$c} eq $f[$i]){
      $sum+=$f[$i+2]*$f[$i+1];
      $weight_sum+=$f[$i+1];
      if($f[$i+3] eq "+"){
         $fwd_sum+=$f[$i+1];
      }else{
         $rev_sum+=$f[$i+1];
      }
    }
  }
  my $mean_pos=$sum/$weight_sum;
  $sum=0;
  #we compute error in position
  for($i=0;$i<$#f;$i+=3){
    if($chrom{$c} eq $f[$i]){
      $sum+=($f[$i+2]-$mean_pos)*($f[$i+2]-$mean_pos)*$f[$i+1];
    }
  }
  $start_pos{$c}=int($mean_pos-$ctg_len{$c}/2);
  $end_pos{$c}=$start_pos{$c}+$ctg_len{$c};
  $error_pos{$c}=$sum/$weight_sum;
  $error_interval{"$start_pos{$c} $end_pos{$c}"}=$error_pos{$c};
  $percent=int($ctg_len{$c}/$ref_chr_len{$chrom{$c}}*100);
  if($fwd_sum>=$rev_sum){
    print "$start_pos{$c} $end_pos{$c} | 1 $ctg_len{$c} | $ctg_len{$c} $ctg_len{$c} | 100 | $ref_chr_len{$chrom{$c}} $ctg_len{$c} | $percent 100 | $chrom{$c} $c\n";
  }else{
    print "$start_pos{$c} $end_pos{$c} | $ctg_len{$c} 1 | $ctg_len{$c} $ctg_len{$c} | 100 | $ref_chr_len{$chrom{$c}} $ctg_len{$c} | $percent 100 | $chrom{$c} $c\n";
  }
}

