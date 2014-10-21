#!/usr/bin/env perl
#
#this code produces joined mega-reads
#
#first we read in PB sequences
my $pbseqfile=$ARGV[0];
my $max_gap=$ARGV[1];
my $allowed_gaps=$ARGV[2];
my $kmer=$ARGV[3];
my $bad_pb=$ARGV[4];
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
    $bad_pb{$pb}="$badstart $badend";
}

my @lines=();
my $outread="";
#now we process the pb+mega-reads file
while($line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
	if(@lines){
	    @lines_sorted = sort by_first_number @lines;
	    $outread=&process_sorted_lines(@lines_sorted);
	    if(not($outread eq "")){
		$indx=0;
		@f=split(/N/,$outread);
		for($i=0;$i<=$#f;$i++){
		    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=400);
		    $indx++;
		}
	    }
	@lines=();
	}
	($rn,$junk)=split(/\s+/,substr($line,1));
    }else{
	push(@lines, $line);
    }
}
#do not forget the last one
        if(@lines){
            @lines_sorted = sort by_first_number @lines;
            $outread=&process_sorted_lines(@lines_sorted);
            if(not($outread eq "")){
                $indx=0;
                @f=split(/N/,$outread);
                for($i=0;$i<=$#f;$i++){
                    print ">$rn.${indx}_",length($f[$i]),"\n$f[$i]\n" if(length($f[$i])>=400);
                    $indx++;
                }
            }
        }



sub by_first_number{
($fn1,$rest)=split(/\s+/,$a);
($fn2,$rest)=split(/\s+/,$b);
return($a <=> $b);
}

sub process_sorted_lines{
    my $outread="";
    $last_coord =-1000000000;
    foreach $l(@_){
        ($bgn,$end,$mbgn,$mend,$mlen,$pb,$mseq,$name)=split(/\s+/,$l);
	$seq=substr($mseq,$mbgn,$mend-$mbgn+1);
        die("inconsistent sequence length") if(not(length($mseq)==$mlen));
        die("pacbio read $pb does not exist in the sequence file!!!") if(not(defined($pbseq{$pb})));

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
                my $max_gap_local;

		if(defined($bad_pb{$pb})){
                my ($bad_start,$bad_end)=split(/\s+/,$bad_pb{$pb});
                $join_allowed=0 if($last_coord<=$bad_start && $bad_start<=$bgn);
                $join_allowed=0 if($last_coord<=$bad_end && $bad_end<=$bgn);
                $join_allowed=0 if($last_coord>=$bad_start && $bgn<=$bad_end);
                }

		$max_gap_local=(length($outread)>length($seq)?length($outread):length($seq));

                if($bgn-$last_coord<$max_gap_local && $join_allowed){#then put N's and later split
		    $outread.=lc(substr($pbseq{$pb},$last_coord+1,$bgn-$last_coord)).$seq;
                }else{
		    $outread.="N".$seq;
                }
            }else{#overlapping
		$join_allowed=1 if($last_mr eq $name); #allow rejoining broken megareads
		#here we check for overlaps
	        my $offset; 
		if($last_coord-$bgn > 10){
                    $ind=index($outread,substr($seq,0,31),length($outread)-($last_coord-$bgn)*1.2);
                    if($ind==-1 || abs(($last_coord-$bgn)-(length($outread)-$ind))>(0.2*($last_coord-$bgn)+10)){
                        $offset=$last_coord-$bgn+1;
                        if($offset > $kmer){
				$join_allowed=0;
			}else{
				$join_allowed=1;
			}
		    }else{
			$offset=length($outread)-$ind;
			$join_allowed=1;
		    }
		}

		if(defined($bad_pb{$pb})){
                my ($bad_start,$bad_end)=split(/\s+/,$bad_pb{$pb});
                $join_allowed=0 if($last_coord>=$bad_start && $bad_start>=$bgn);
                $join_allowed=0 if($last_coord>=$bad_end && $bad_end>=$bgn);
                $join_allowed=0 if($bgn>=$bad_start && $last_coord<=$bad_end);
                }

#print "$join_allowed $last_coord $bgn INDEX $ind ",length($outread)," ",$last_coord-$bgn+1," ",length($outread)-$ind,"\n";
		if($join_allowed){
		    $outread.=substr($seq,$offset);
		}else{
		    $outread.="N".$seq;
		}
            }
        }
	$last_coord=$end;
	$last_mr=$name;
	$last_len=length($seq);
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

