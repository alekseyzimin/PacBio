#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
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
    my $min_match=35;

    foreach $l(@args){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=@{$l};
        if($bgn<=$last_coord && $last_coord-$bgn<=$min_match){#if small overlap, we extend the matches to find a better overlap
	    my $tlen=length($last_tail);
	    if($tlen<$min_match){
		$outseq.=$last_tail;
		$last_coord+=$tlen;
	    }else{
		$outseq.=substr($last_tail,0,$min_match);
		$last_coord+=$min_match;
	    }
	    if($mbgn<$min_match){
		$mbgn=1;
		$bgn-=$mbgn;
		if($bgn<1){
		    $mbgn-=($bgn-1);
		    $bgn=1;
		}
	    }
        }
	$seq=substr($mseq,$mbgn-1,$mend-$mbgn+1);
        die("inconsistent sequence length") if(not(length($mseq)==$mlen));
        die("pacbio read $pb does not exist in the sequence file!!!") if(not(defined($pbseq{$pb})));
        $gap_index++;
        if($outread eq ""){
	    $outread=$seq; # the first chunk
        }else{
            my @k1s=split(/_/,$last_mr);
            my @k2s=split(/_/,$name);
            $k1s[$#k1s] =substr($k1s[$#k1s],0,-1);
            $k2s[0] = substr($k2s[0],0,-1);
            $str="$pb $k1s[$#k1s] $k2s[0]";
            $str="$pb $k2s[0] $k1s[$#k1s]" if($k1s[$#k1s]>$k2s[0]);

            $join_allowed=0;
            if($bgn<=$last_coord){#overlapping
		# we now allowing this globally $join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads
		my $ind=-1;
		my %ind=();
                my $offset=-1;

		if($last_coord-$bgn >= $min_match){ #it is possible to check for overlap
		    for(my $j=0;$j<10;$j++){
                    	my $ttt=index($outread,substr($seq,$j,$min_match),length($outread)-($last_coord-$bgn)*$fudge_factor);
		        $ind{$ttt-$j}++ if($ttt>-1);
		    }
		    my $max_ind=-1;
		    foreach my $ttt (keys %ind){
			if($ind{$ttt}>$max_ind){$max_ind=$ind{$ttt};$ind=$ttt;}
		    }
                    if($ind==-1 || ($ind>-1 && abs(($last_coord-$bgn)-(length($outread)-$ind))>(0.3*($last_coord-$bgn)+10))){ #if no overlap or overlap inconsistent with implied
		        $join_allowed=0;
                    }else{
                        $join_allowed=1;
                    }
		}

		if($join_allowed){#here if allowed means that either the overlap was too short or match was found
		    if($offset>-1){
			$outread=substr($outread,0,length($outread)-$offset).$seq;
		    }elsif($ind>-1){
			$outread=substr($outread,0,$ind).$seq;
		    }else{
                        #we should never get here
                        die("error in joining $offset $ind $pb $name");
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
	$last_tail=(length($mseq)>$mend) ? "" : substr($mseq,$mend+1);
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

