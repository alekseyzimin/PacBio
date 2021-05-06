#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use mummer;

my $allowed_gaps=$ARGV[0];
my $max_gap=$ARGV[1];
my $min_len_output=500;
my $fudge_factor=1.2;
my $pbseq="";

open(FILE,$allowed_gaps);
while($line=<FILE>){
    chomp($line);
    @f=split(/\s+/,$line);
    $allowed{"$f[0] $f[2] $f[3]"}=$f[-1] if(defined($f[3]));
}

my @lines=();
my $outread="";
#now we process the pb+mega-reads file
while($line=<STDIN>){
  chomp($line);
  if($line =~ /^>/){
    if(@lines && not($rn eq "") && not($pbseq eq "")){
      $outread = "";
      $outread = process_sorted_lines(sort {$$a[0] <=> $$b[0]} @lines);
      if(not($outread eq "")){
        $indx=0;
        @f=split(/(N{1,})/,$outread);
        if(scalar(@f)==1){
          print ">$rn.1_",length($outread),"\n$outread\n" if(length($outread)>=$min_len_output);
        }else{
          for($i=0;$i<=$#f;$i+=2){
            print STDERR ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=$min_len_output);
            $indx+=length($f[$i]);
            $indx+=length($f[$i+1]) if($f[$i]<$#f);
          }
        }
      }
      @lines=();
    }
    ($rn,$pbseq)=split(/\s+/,substr($line,1));
  }else{
    my @ttt=split(/\s+/,$line);
    push(@lines, \@ttt) if($#ttt==7);
  }
}
#do not forget the last one
if(@lines && not($rn eq "") && not($pbseq eq "")){
  $outread = process_sorted_lines(sort {$$a[0] <=> $$b[0]} @lines);
  if(not($outread eq "")){
    $indx=0;
    @f=split(/(N{1,})/,$outread);
    if(scalar(@f)==1){
      print ">$rn.1_",length($outread),"\n$outread\n" if(length($outread)>=$min_len_output);
    }else{
      for($i=0;$i<=$#f;$i+=2){
        print STDERR ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=$min_len_output);
        $indx+=length($f[$i]);
        $indx+=length($f[$i+1]) if($f[$i]<$#f);
      }
    }
  }
}


sub process_sorted_lines{
    my $outread="";
    $last_coord =-1000000000;
    $last_mend=0;
    my @max_gap_local_fwd=();
    my @max_gap_local_rev=();
    my @args=@_;
    my $gap_coeff=1;
    my $outread_len=0;
    my $seq_len=0;
    my $sum_chunk_size=0;
    my $num_chunks=0;
    my $min_match=17;

    for(my $i=0;$i<=$#args;$i++){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
	$sum_chunk_size+=($mend-$mbgn);
	$num_chunks++;
    }

    return($outread) if($sum_chunk_size/$num_chunks<$min_len_output);#average chunk size must be >$min_len_output

    for(my $i=0;$i<=$#args;$i++){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
	$outread_len+=$mend-$mbgn;
	$seq_len=$mend-$mbgn;
        $max_gap_local=$gap_coeff*($outread_len>$seq_len?$outread_len:$seq_len);
        $max_gap_local=$max_gap if($max_gap_local>$max_gap);
	push(@max_gap_local_fwd,$max_gap_local);
    }

    $outread_len=0;
    $seq_len=0;
    for(my $i=$#args;$i>=0;$i--){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
        $outread_len+=$mend-$mbgn;
        $seq_len=$mend-$mbgn;
        $max_gap_local=$gap_coeff*($outread_len>$seq_len?$outread_len:$seq_len);
        $max_gap_local=$max_gap if($max_gap_local>$max_gap);
        unshift(@max_gap_local_rev,$max_gap_local);
    }
    
    my $gap_index=-1;
    foreach $l(@args){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$l};
	$seq=substr($mseq,$mbgn-1,$mend-$mbgn+1);
        next if(not(length($mseq)==$mlen));
        #die("inconsistent sequence length for $pb") if(not(length($mseq)==$mlen));
        die("long read $pb does not exist in the sequence file!!!") if(not(defined($pbseq)));
        $gap_index++;
        if($outread eq ""){
	    $outread=$seq; # the first chunk
        }else{
            next if($end<=$last_coord);#we should not allow containment
            my @k1s=split(/_/,$last_mr);
            my @k2s=split(/_/,$name);
            $k1s[$#k1s] =substr($k1s[$#k1s],0,-1);
            $k2s[0] = substr($k2s[0],0,-1);
            $str="$pb $k1s[$#k1s] $k2s[0]";
            $str="$pb $k2s[0] $k1s[$#k1s]" if($k1s[$#k1s]>$k2s[0]);

            $join_allowed=0;
	    $join_allowed=$allowed{$str} if(defined($allowed{$str})); #allow joins that are in multiple pacbios with a gap.  without gap we need to see at least 11 bp overlap/alignment 
            $join_allowed=1 if($last_mr eq $name && $bgn-$last_coord<-5); #allow rejoining broken megareads when overlapping ends

            if($bgn>$last_coord){#if gap -- check if the closure is allowed
		#$max_gap_local=$max_gap_local_fwd[$gap_index]<$max_gap_local_rev[$gap_index]?$max_gap_local_fwd[$gap_index]:$max_gap_local_rev[$gap_index];
                $max_gap_local=$max_gap;
                $max_gap_local=$max_gap_local/2 if($join_allowed==-1); 
                if($bgn-$last_coord<=$max_gap_local && ($join_allowed==1 || $join_allowed==-1)){#join or then put N's and later split
		    $outread.=lc(substr($pbseq,$last_coord,$bgn-$last_coord-1)).$seq;
                }else{
		    $outread.="N"x($bgn-$last_coord).$seq;
                }
            }else{#overlapping
# we now allowing this globally $join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads
              my $ind=-1;
              $join_allowed=abs($join_allowed);
              my $slack=int(($last_coord-$bgn)*0.05)+$min_match;
              my @args=@_;
              my $o=mummer::Options->new;
              $o->minmatch($min_match);
              $o->mincluster($min_match);
              $o->forward();
              my $overlap=$last_coord-$bgn+$slack;
              my $ind2=length($outread)-$overlap+$slack-1;#default implied overlap
              #print "DEBUG: overlap $overlap slack $slack join_allowed $join_allowed\n";
              #print "DEBUG:\n>o\n",substr($outread,length($outread)-$overlap),"\n>s\n",substr($seq,0,$overlap),"\n";
              my $a = mummer::align_sequences(substr($outread,length($outread)-$overlap),substr($seq,0,$overlap), $o);
              if(scalar(@$a)>0){
                #print $$a[0]{sA}," ",$$a[0]{eA}," ",$$a[0]{sB}," ",$$a[0]{eB}," $$a[0]{Errors} $$a[0]{SimErrors} $$a[0]{NonAlphas}\n";
                $seq=substr($seq,$$a[0]{sB}-1);
                $ind=length($outread)-$overlap+$$a[0]{sA}-1;
              }
              #print "DEBUG ind $ind ind2 $ind2\n";
              if($ind==-1){
                if($overlap<50 && $join_allowed){#small overlap, probably repeat, nucmer did not find it, join is allowed
                  $ind=$ind2;
                  $join_allowed=1;
                }else{
                  $join_allowed=0;
                }
              }else{
                $join_allowed=1;
              }

              if($join_allowed==1){#here if allowed means that either the overlap was too short or match was found
                if($ind>-1){
                  $outread=substr($outread,0,$ind).$seq;
                  #print "DEBUG: $outread\n";
                }else{
#we should never get here
                  die("error in joining $offset $ind $pb $name $join_allowed");
                }
              }else{
                $outread.="N".$seq;
              }
            }
        }
        $last_coord=$end;
        $last_mend=$mend;
        $last_mr=$name;
    }
    return($outread);
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

