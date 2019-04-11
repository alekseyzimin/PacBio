#!/usr/bin/env perl
#first we read in sequences of contigs to merge
my $contigs_to_merge=$ARGV[0];
my $ctg="",$seq="";
my %seq =(),%output=(),@ctgnames=();
open(FILE,$contigs_to_merge);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $seq{$ctg}=$seq;
#to preserve order
      push(@ctgnames,$ctg);
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $seq{$ctg}=$seq;
  push(@ctgnames,$ctg);
}

#then we read in the merging sequences and trim points; we record the sequence that follows each contig
my $merges=$ARGV[1];
open(FILE,$merges);
while($line=<FILE>){
    chomp($line);
    my ($c1,$t1,$d1,$c2,$t2,$d2,$g,$s)=split(/\s+/,$line);

    if($d1 eq "F"){
      $trim3{$c1}=$t1;
    }else{
      $trim5{$c1}=$t1;
    }

    if($d2 eq "F"){
      $trim5{$c2}=$t1;
    }else{
      $trim3{$c1}=$t1;
    }   

    $gseq{"$c1$d1$c2$d2"}=$s;
    $d1 =~ tr/FR/RF/;
    $d2 =~ tr/FR/RF/;
    $gseq{"$c2$d2$c1$d1"}=reverse_complement($s);
} 

#trim chunks
foreach my $ctg(keys %seq){
my $beg=0;
my $end=length($seq{$ctg});
$beg+=$trim5{$ctg} if(defined($trim5{$ctg}));
$end-=$trim3{$ctg} if(defined($trim3{$ctg}));
die("trim $beg for contig $ctg outside range") if($beg>= length($seq{$ctg}));
die("trim $end for contig $ctg outside range") if($end<=0);
$seq{$ctg}=substr($seq{$ctg},$beg,$end-$beg);
}

#now read in the merges file
while($line=<STDIN>){
    chomp($line);
    my @f=split(/\s+/,$line);
    my $readname=$f[0];
    for(my $i=1;$i<=$#f;$i++){
      if($i%3==0){
        my @ff=split(/\./,$f[$i]);
        $readname.=":$ff[-1]";
      }else{
        $readname.=":$f[$i]";
      }
    }
    print ">$readname\n";
    my $merge_arg="$f[-3]$f[-2]$f[$i]$f[$i+1]";
#output first chunk
    if($f[1] eq "R"){
	print reverse_complement($seq{$f[0]});
    }else{
	print $seq{$f[0]};
    }
    $output{$f[0]}=1;
    for($i=3;$i<$#f;$i+=3){
        $merge_arg="$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]";
	if($f[$i-1]>0){
            die("gap $merge_arg not found") if(not(defined($gseq{$merge_arg})));
            die("sequence $f[$i] not found") if(not(defined($seq{$f[$i]})));
            print $gseq{$merge_arg};
	    if($f[$i+1] eq "R"){
		print reverse_complement($seq{$f[$i]});
	    }else{
		print $seq{$f[$i]};
	    }
	}else{#negative gap
	    if($f[$i+1] eq "R"){
		print substr(reverse_complement($seq{$f[$i]}),-$f[$i-1]);
	    }else{
		print substr($seq{$f[$i]},-$f[$i-1]);
	    }
	}
	$output{$f[$i]}=1;
    }
    print "\n";
}

foreach $c(@ctgnames){
    print ">$c\n$seq{$c}\n" if(not(defined($output{$c})));
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

