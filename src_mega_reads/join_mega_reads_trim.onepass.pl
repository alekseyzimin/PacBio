#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $pbseqfile=$ARGV[0];
my $allowed_gaps=$ARGV[1];
my $kmer=$ARGV[2];
my $bad_pb=$ARGV[3];
my $fudge_factor=1.2;
my $min_len_output=400;
$kmer*=$fudge_factor;

my $rn="";
my %pbseq;
open(FILE,$pbseqfile);
while($line=<FILE>){
    chomp($line);
    if($line =~ /^>/){
	@f=split(/\s+/,$line);
	$rn=substr($f[0],1);
    }else{
	$pbseq{$rn}.=$line;
    }
}

open(FILE,$allowed_gaps);
while($line=<FILE>){
    chomp($line);
    @f=split(/\s+/,$line);
    $allowed{"$f[0] $f[2] $f[3]"}=$f[-1];
}

open(FILE,$bad_pb);
while($line=<FILE>){
    chomp($line);
    ($pb,$badstart,$badend)=split(/\s+/,$line);
    $bad_pb{$pb}.="$badstart $badend ";
}

my @lines=();
my $outread="";
#now we process the pb+mega-reads file
while($line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
	if(@lines){
	    $outread = "";
	    $outread = process_sorted_lines(sort {$$a[0] <=> $$b[0]} @lines) if($#lines<100);#no more than 100 chunks per PB read
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

    for(my $i=0;$i<=$#args;$i++){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$args[$i]};
	$sum_chunk_size+=($mend-$mbgn);
	$num_chunks++;
	}

    return($outread) if($sum_chunk_size/$num_chunks<500);#average chunk size must be >500bp

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
	$ibgn=$bgn-$mbgn+1;
	$iend=$end+($mlen-$mend);
	$seq=$mseq;
        die("inconsistent sequence length") if(not(length($mseq)==$mlen));
        die("pacbio read $pb does not exist in the sequence file!!!") if(not(defined($pbseq{$pb})));
        $gap_index++;
        if($outread eq ""){
	    $outread=substr($seq,$mbgn-1);#the first one is trimmed at beg, but allowed full length
        }else{
            $join_allowed=0;

            my @k1s=split(/_/,$last_mr);
            my @k2s=split(/_/,$name);
            $k1s[$#k1s] =substr($k1s[$#k1s],0,-1);
            $k2s[0] = substr($k2s[0],0,-1);
            $str="$pb $k1s[$#k1s] $k2s[0]";
            $str="$pb $k2s[0] $k1s[$#k1s]" if($k1s[$#k1s]>$k2s[0]);

	    if(defined($allowed{$str})){
		$join_allowed=$allowed{$str};
	    }else{
		$join_allowed=1;
	    }

            if(defined($bad_pb{$pb})){
                    my @bad_coords=split(/\s+/,$bad_pb{$pb});
                    for($i=0;$i<=$#bad_coords;$i+=2){
                        my $bad_start=$bad_coords[$i];
                        my $bad_end=$bad_coords[$i+1];
                        $join_allowed=0 if($last_coord<=$bad_start && $bad_start<=$bgn);
                        $join_allowed=0 if($last_coord<=$bad_end && $bad_end<=$bgn);
                        $join_allowed=0 if($last_coord>=$bad_start && $bgn<=$bad_end);
                    }
                }


            if($ibgn>$last_icoord){#if implied gap -- check if the closure is allowed
		$max_gap_local=$max_gap_local_fwd[$gap_index]<$max_gap_local_rev[$gap_index]?$max_gap_local_fwd[$gap_index]:$max_gap_local_rev[$gap_index];
                if($ibgn-$last_icoord<$max_gap_local && $join_allowed){#then allow raw sequence
                    my $fill_seq=lc(substr($pbseq{$pb},$last_icoord,$ibgn-$last_icoord-1));
		    $outread.=$fill_seq.$seq;
                    if($k1s[$#k1s]>$k2s[0]){
			print STDERR "$pb $k2s[0] $k1s[$#k1s] ",$last_icoord," real ",$bgn-$last_coord-1," impl ",$ibgn-$last_icoord-1," ",reverse_complement($fill_seq),"\n";
			}else{
			print STDERR "$pb $k1s[$#k1s] $k2s[0] ",$last_icoord," real ",$bgn-$last_coord-1," impl ",$ibgn-$last_icoord-1," $fill_seq\n";
			}		
                }else{
		    $outread=substr($outread,0,length($outread)-($last_icoord-$last_coord))."N"x($ibgn-$last_icoord).substr($seq,$mbgn-1);
                }
            }elsif($ibgn<=$last_icoord && $bgn>$last_coord){#implied overlap but no alignment overlap
                if($bgn-$last_coord<$max_gap_local && $join_allowed){#then allow raw sequence
                    my $fill_seq=lc(substr($pbseq{$pb},$last_coord,$bgn-$last_coord-1));
                    $outread=substr($outread,0,length($outread)-($last_icoord-$last_coord)).$fill_seq.substr($seq,$mbgn-1);
                    if($k1s[$#k1s]>$k2s[0]){
                        print STDERR "$pb $k2s[0] $k1s[$#k1s] ",$last_coord," real ",$bgn-$last_coord-1," impl ",$ibgn-$last_icoord-1," ",reverse_complement($fill_seq),"\n";
                        }else{
                        print STDERR "$pb $k1s[$#k1s] $k2s[0] ",$last_coord," real ",$bgn-$last_coord-1," impl ",$ibgn-$last_icoord-1,"  $fill_seq\n";
                        }
                }else{
                    $outread=substr($outread,0,length($outread)-($last_icoord-$last_coord))."N"x($bgn-$last_coord).substr($seq,$mbgn-1);
                }
		}else{#we have a real overlap -- forget about the implied one,here we must first trim the outread and the sequence
		$outread=substr($outread,0,length($outread)-($last_icoord-$last_coord));
                $seq=substr($seq,$mbgn-1);
                $join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads
                my $min_match=25;
                my $ind=-1;
                my %ind=();

                if($last_coord-$bgn > $min_match){
                    for(my $j=0;$j<10;$j++){
                        my $ttt=index($outread,substr($seq,$j,$min_match),length($outread)-($last_coord-$bgn)*$fudge_factor);
                        $ind{$ttt-$j}++ if($ttt>-1);
                        }
                    my $max_ind=-1;
                    foreach my $ttt (keys %ind){
                        if($ind{$ttt}>$max_ind){$max_ind=$ind{$ttt};$ind=$ttt;}
                        }

                    if($ind==-1 || ($ind>-1 && abs(($last_coord-$bgn)-(length($outread)-$ind))>(0.2*($last_coord-$bgn)+10))){
                        $offset=$last_coord-$bgn+1;
                        if($offset > $kmer){
                            $join_allowed=0;
                        }else{
                            $join_allowed=1;
                        }
                    }else{
                        $join_allowed=1;
                    }
                }

                if($join_allowed){
                    if($ind==-1){
                    $outread=substr($outread,0,length($outread)-($last_coord-$bgn+1)).$seq;
                        }else{
                    $outread=substr($outread,0,$ind).$seq;
                        }
                }else{
                    $outread.="N".$seq;
                }
            }
        }
	$last_coord=$end;
	$last_icoord=$iend;
	$last_mr=$name;
	$last_len=$mlen;
        last if($last_coord>=length($pbseq{$pb}));	
    }
    return(substr($outread,0,length($outread)-($last_icoord-$last_coord)));
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

