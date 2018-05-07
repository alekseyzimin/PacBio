#!/usr/bin/env perl
#this code simply merges contig sequences based on output of show-coords -o
#each merge line must be formatted like this
#       1     8287  |    55192    63478  |     8287     8287  |   100.00  |  1871264    63478  |  1  1  jcf7190000000040 jcf7190000012554        [BEGIN]
#

$seq_file=$ARGV[0];

open(FILE,$seq_file);
while($line=<FILE>){
  chomp($line);
  if($line =~ /^>/){
    @f=split(/\s+/,$line);
    $sn=substr($f[0],1);
  }else{
    $seq{$sn}.=$line;
  }
}

$slack=500;
while($line=<STDIN>){
chomp($line);
$line=~s/^\s+//;
@f=split(/\s+/,$line);
if(defined($seq{$f[16]}) && defined($seq{$f[17]})){
  if($f[0]<$slack){
#merge off the begin of the ref contig 
    if($f[3]<$f[4]){
      if($f[4]>$f[12]-$slack){
        $new_contig=$f[17]."_F_".$f[16]."_F";
        $new_seq=substr($seq{$f[17]},0,$f[3]-1).substr($seq{$f[16]},$f[0]-1);
        delete($seq{$f[17]});
        delete($seq{$f[16]});
        $seq{$new_contig}=$new_seq;
      }else{
        print STDERR "WARNING ahang too big in $line\n"; 
      }
    }else{
      if($f[4]<$slack){
        $new_contig=$f[17]."_R_".$f[16]."_F";
        $new_seq=reverse_complement(substr($seq{$f[17]},$f[3])).substr($seq{$f[16]},$f[0]-1);
        delete($seq{$f[17]});
        delete($seq{$f[16]});
        $seq{$new_contig}=$new_seq;
      }else{
        print STDERR "WARNING bhang too big in $line\n";
      }
    }
  }elsif($f[1]>$f[11]-$slack){
#merge off the end of the ref contig
    if($f[3]<$f[4]){
      if($f[3]<$slack){
        $new_contig=$f[16]."_F_".$f[17]."_F";
        $new_seq=substr($seq{$f[16]},0,$f[0]).substr($seq{$f[17]},$f[3]);
        delete($seq{$f[17]});
        delete($seq{$f[16]});
        $seq{$new_contig}=$new_seq;
      }else{
        print STDERR "WARNING ahang too big in $line\n";
      }
    }else{ 
      if($f[3]>$f[12]-$slack){
        $new_contig=$f[16]."_F_".$f[17]."_R";
        $new_seq=substr($seq{$f[16]},0,$f[0]).reverse_complement(substr($seq{$f[17]},0,$f[3]-1));
        delete($seq{$f[17]});
        delete($seq{$f[16]});
        $seq{$new_contig}=$new_seq;
      }else{
        print STDERR "WARNING bhang too big in $line\n";
      }
    }
  }else{
    print STDERR "WARNING ahang too big in $line\n";
  }
}else{
  print STDERR "WARNING contigs $f[16] and $f[17] already merged!\n";

}
}

foreach $s(keys %seq){
print ">$s\n$seq{$s}\n";
}


sub reverse_complement{
  my $str=$_[0];
  $str =~ tr/acgtACGTNn/tgcaTGCANn/;
  $str = reverse ($str);
  return ($str);
}

