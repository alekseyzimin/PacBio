#!/usr/bin/env perl
#first we read in sequences of contigs to merge
my $contigs_to_merge=$ARGV[0];
my $ctg="",$seq="";
my %seq =(),%output=();
open(FILE,$contigs_to_merge);
while($line=<FILE>){
    chomp($line);
    if($line=~/^\>/){
	my @f=split(/\s+/,$line);
	if(not($seq eq "")){
	    $seq{$ctg}=$seq;
	}
	$ctg=substr($f[0],1);
	$seq="";
    }else{
	$seq.=$line;
    }
}
if(not($seq eq "")){
    $seq{$ctg}=$seq;
}

#then we read in the merging sequences; we record the sequence that follows each contig
my $merges=$ARGV[1];
open(FILE,$merges);
while($line=<FILE>){
    #print STDERR $line;
    chomp($line);
    my ($c1,$o1,$d1,$c2,$o2,$d2,$g,$s)=split(/\s+/,$line);
    $gseq{"$c1$d1$c2$d2"}=$s;
    $oh1{"$c1$d1$c2$d2"}=$o1;
    $oh2{"$c1$d1$c2$d2"}=$o2;
    $d1 =~ tr/FR/RF/;
    $d2 =~ tr/FR/RF/;
    $gseq{"$c2$d2$c1$d1"}=reverse_complement($s);
    $oh1{"$c2$d2$c1$d1"}=$o2;
    $oh2{"$c2$d2$c1$d1"}=$o1;
} 


#now read in the merges file
while($line=<STDIN>){
    chomp($line);
    my $outlen=0,$sumlen=0;
    #print "DEBUG $line\n";
    my @f=split(/\s+/,$line);
    print ">",join(":",@f),"\n";
#output first contig
    $oh2=$oh1{"$f[0]$f[1]$f[3]$f[4]"};
    $len=length($seq{$f[0]});
    $sumlen+=$len;
    #print "DEBUG $f[0] $len $oh2\n";
    if($f[1] eq "R"){
	print substr(reverse_complement($seq{$f[0]}),0,$len-$oh2);
        #print "\nDEBUG ",length(substr(reverse_complement($seq{$f[0]}),0,$len-$oh2)),"\n";
        $outlen+=length(substr(reverse_complement($seq{$f[0]}),0,$len-$oh2));
    }else{
	print substr($seq{$f[0]},0,$len-$oh2);
        #print "\nDEBUG ",length(substr($seq{$f[0]},0,$len-$oh2)),"\n";
        $outlen+=length(substr($seq{$f[0]},0,$len-$oh2));
    }
    $output{$f[0]}=1;
#now the same for the rest we first output the previous gap (or trim) and then the contig
    for($i=3;$i<$#f;$i+=3){
      $oh1=$oh2{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"};
      $oh1=0 if($oh1<0);
      $oh2=0; 
      $oh2=$oh1{"$f[$i]$f[$i+1]$f[$i+3]$f[$i+4]"} if($i+4<=$#f);
      $len=length($seq{$f[$i]});
      $sumlen+=$len;
      #print "\nDEBUG $i $f[$i] ",$f[$i-1]," $len $oh1 $oh2\n";
	if($f[$i-1]>0){
            die("gap $f[$i-3]$f[$i-2]$f[$i]$f[$i+1] not found") if(not(defined($gseq{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"})));
            die("sequence $f[$i] not found") if(not(defined($seq{$f[$i]})));
            print $gseq{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"};
            #print "\nDEBUG gap len ",length($gseq{"$f[$i-3]$f[$i-2]$f[$i]$f[$i+1]"})," $oh1 $f[$i] $oh2 $f[$i+1]\n";
	    if($f[$i+1] eq "R"){
		print substr(reverse_complement($seq{$f[$i]}),$oh1,$len-$oh1-$oh2);
                #print "\nDEBUG ",length(substr(reverse_complement($seq{$f[$i]}),$oh1,$len-$oh1-$oh2)),"\n";
                $outlen+=length(substr(reverse_complement($seq{$f[$i]}),$oh1,$len-$oh1-$oh2));
	    }else{
		print substr($seq{$f[$i]},$oh1,$len-$oh1-$oh2);
                #print "\nDEBUG ",length(substr($seq{$f[$i]},$oh1,$len-$oh1-$oh2)),"\n";
                $outlen+=length(substr($seq{$f[$i]},$oh1,$len-$oh1-$oh2));
	    }
	}else{#negative gap
            #print "\nDEBUG negative gap len ",$f[$i-1]," $oh1 $f[$i] $oh2 $f[$i+1]\n";
	    if($f[$i+1] eq "R"){
		print substr(substr(reverse_complement($seq{$f[$i]}),$oh1,$len-$oh1-$oh2),-$f[$i-1]);
                #print "\nDEBUG ",length(substr(substr(reverse_complement($seq{$f[$i]}),$oh1,$len-$oh1-$oh2),-$f[$i-1])),"\n";
                $outlen+=length(substr(substr(reverse_complement($seq{$f[$i]}),$oh1,$len-$oh1-$oh2),-$f[$i-1]));
	    }else{
		print substr(substr($seq{$f[$i]},$oh1,$len-$oh1-$oh2),-$f[$i-1]);
                #print "\nDEBUG ",length(substr(substr($seq{$f[$i]},$oh1,$len-$oh1-$oh2),-$f[$i-1])),"\n";
                $outlen+=length(substr(substr($seq{$f[$i]},$oh1,$len-$oh1-$oh2),-$f[$i-1]));
	    }
	}
	$output{$f[$i]}=1;
    }
    #print "DEBUG output $outlen $sumlen\n";
    print "\n";
}

foreach $c(keys %seq){
   print ">$c\n$seq{$c}\n" if(not(defined($output{$c})));
}

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}

