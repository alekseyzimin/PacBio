#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $kmer=40;
my $fudge_factor=1.2;
my $min_len_output=400;
$kmer*=$fudge_factor;

my @lines=();
my $outread="";
#now we process the pb+mega-reads file
while($line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
	if(@lines){
	    $outread = "";
	    $outread = process_sorted_lines(sort {$$a[0] <=> $$b[0]} @lines);
	    if(not($outread eq "")){
		$indx=0;
		@f=split(/(N{1,})/,$outread);
		for($i=0;$i<=$#f;$i+=2){
		    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=$min_len_output);
		    $indx+=length($f[$i]);
		    $indx+=length($f[$i+1]) if($f[$i]<$#f);
		}
	    }
	    @lines=();
	}
	($rn,$junk)=split(/\s+/,substr($line,1));
    }else{
	my @ttt=split(/\s+/,$line);
	push(@lines, \@ttt);
    }
}
#do not forget the last one
if(@lines){
    $outread = process_sorted_lines(sort {$$a[0] <=> $$b[0]} @lines);
    if(not($outread eq "")){
	$indx=0;
	@f=split(/(N{1,})/,$outread);
	for($i=0;$i<=$#f;$i+=2){
	    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=$min_len_output);
	    $indx+=length($f[$i]);
	    $indx+=length($f[$i+1]) if($f[$i]<$#f);
	}
    }
}


sub process_sorted_lines{
    my $outread="";
    my $last_seq="";
    $last_coord =-1000000000;
    my @max_gap_local_fwd=();
    my @max_gap_local_rev=();
    my @args=@_;
    my $max_gap=750;
    my $gap_coeff=1;
    my $outread_len=0;
    my $seq_len=0;
    my $sum_chunk_size=0;
    my $num_chunks=0;
    my $join_allowed=0;

    for(my $i=0;$i<=$#args;$i++){
      ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
      $sum_chunk_size+=($mend-$mbgn);
      $num_chunks++;
    }

    return($outread) if($sum_chunk_size/$num_chunks<400);#average chunk size must be >400bp

    my $gap_index=-1;
    foreach $l(@args){
      ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$l};
      $seq=substr($mseq,$mbgn-1,$mend-$mbgn+1);
      die("inconsistent sequence length") if(not(length($mseq)==$mlen));
      $gap_index++;
      if($outread eq ""){
        $outread=$seq; # the first chunk
      }else{
        next if($end<=$last_coord);
        $join_allowed=1 if($last_mr eq $name && $bgn-$last_coord<5); #allow rejoining broken megareads when overlapping ends/gap small


# we now allowing this globally $join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads
        my $min_match=25;
        my $ind=-1;
        my %ind=();
        my $offset=-1;

        if($last_coord-$bgn > $min_match){ #it is possible to check for overlap
          for(my $j=0;$j<10;$j++){
            my $ttt=index($outread,substr($seq,$j,$min_match),length($outread)-($last_coord-$bgn)*$fudge_factor);
            $ind{$ttt-$j}++ if($ttt>-1);
          }
          my $max_ind=-1;
          foreach my $ttt (keys %ind){
            if($ind{$ttt}>$max_ind){$max_ind=$ind{$ttt};$ind=$ttt;}
          }
          if($ind==-1 || ($ind>-1 && abs(($last_coord-$bgn)-(length($outread)-$ind))>(0.2*($last_coord-$bgn)+10))){ #if no overlap or overlap inconsistent with implied
            $join_allowed=0;
          }else{
            $join_allowed=1;
          }

        if($join_allowed){#here if allowed means that either the overlap was too short or match was found
          if($ind>-1){
            $outread=substr($outread,0,$ind).$seq;
          }else{
#we should never get here
            die("error in joining $offset $ind $pb $name");
          }
        }else{
          $outread.="N".$seq;
        }
      }else{#gap
        $outread.="N".$seq;
      }
      }
      $last_coord=$end;
      $last_mr=$name;
      $last_seq=$seq;
      $last_mend=$mend;
    }
    return($outread);
}

sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}

