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
	    $outread = process_sorted_lines(sort {$$a[0] <=> $$b[0]} @lines);
	    if(not($outread eq "")){
		$indx=0;
		@f=split(/(N{1,})/,$outread);
		for($i=0;$i<=$#f;$i+=2){
		    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=500);
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
                    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=500);
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
    my $max_gap=1250;
    my $gap_coeff=1;
    my $outread_len=0;
    my $seq_len=0;

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
        die("inconsistent sequence length") if(not(length($mseq)==$mlen));
        die("pacbio read $pb does not exist in the sequence file!!!") if(not(defined($pbseq{$pb})));
        $gap_index++;
        if($outread eq ""){
	    $outread=$seq;
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

            if($bgn>$last_coord){#if gap -- check if the closure is allowed

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

		$max_gap_local=$max_gap_local_fwd[$gap_index]<$max_gap_local_rev[$gap_index]?$max_gap_local_fwd[$gap_index]:$max_gap_local_rev[$gap_index];

                if($bgn-$last_coord<$max_gap_local && $join_allowed){#then put N's and later split
		    $outread.=lc(substr($pbseq{$pb},$last_coord,$bgn-$last_coord-1)).$seq;
                }else{
		    $outread.="N"x($bgn-$last_coord).$seq;
                }
            }else{#overlapping
		$join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads
		#here we check for overlap
	 	my $min_match=27;
		my $ind=-1;
		my %ind=();

		if($last_coord-$bgn > $min_match){
		    for(my $j=0;$j<10;$j++){
                    	my $ttt=index($outread,substr($seq,$j,$min_match),length($outread)-($last_coord-$bgn)*1.2);
		        $ind{$ttt-$j}++ if($ttt>-1);
		    	}
		    my $max_ind=-1;
		    foreach my $ttt (keys %ind){
			#print "DEBUG possible ind $ttt freq $ind{$ttt} implied offset ",length($outread)-($last_coord-$bgn),"\n";
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
		}else{
		    $join_allowed=1  unless($allowed{$str}==0);
		}
		
		if(defined($bad_pb{$pb})){
		    my @bad_coords=split(/\s+/,$bad_pb{$pb});
		    for($i=0;$i<=$#bad_coords;$i+=2){
			my $bad_start=$bad_coords[$i];
			my $bad_end=$bad_coords[$i+1];
			$join_allowed=0 if($last_coord>=$bad_start && $bad_start>=$bgn);
			$join_allowed=0 if($last_coord>=$bad_end && $bad_end>=$bgn);
			$join_allowed=0 if($bgn>=$bad_start && $last_coord<=$bad_end);
		    }
		}
#print "$join_allowed $last_coord $bgn INDEX $ind ",length($outread)," ",$last_coord-$bgn+1," ",length($outread)-$ind,"\n";
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
	$last_mr=$name;
	$last_seq=$seq;
	$last_mend=$mend;
        last if($last_coord>=length($pbseq{$pb}));	
    }
    return($outread);
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

