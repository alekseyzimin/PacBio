#!/usr/bin/env perl
#
#this script introduces single base and insertion/deletion errors int fasta file and records them in vcf file
my $ref_contigs=$ARGV[0];
my $error_rate=$ARGV[1];
my $max_indel=20;
$max_indel=$ARGV[2] if($ARGV[2]>0);
my %rseq=(),$seq="",$ctg="";

open(FILE,$ref_contigs);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $rseq{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $rseq{$ctg}=$seq;
}

foreach $ctg (keys %rseq){
$seq=$rseq{$ctg};
for($i=2*$max_indel;$i<length($seq)-2*$max_indel;$i++){
  next if(uc(substr($seq,$i,1)) eq "N");
  if(rand(1)<$error_rate){#add an error
    if(rand(1)<0.9){#error is substitution
      $sub="A";
      $sub="G" if(uc(substr($seq,$i,1)) eq "A");
      print "$ctg\t",$i+1,"\t.\t",substr($seq,$i,1),"\t$sub\t*\t*\t*\t*\t1:1:1:0:0:10:10:0\n";
    }else{#error is an indel
      $size=int(rand($max_indel-1))+1;
      if(rand(1)<0.5){#error is deletion
        print "$ctg\t",$i+1,"\t.\t",substr($seq,$i,$size+1),"\t",substr($seq,$i,1),"\t*\t*\t*\t*\t1:1:1:0:0:10:10:0\n";
      }else{
        #error is insertion
        $sub="A"x$size;
        print "$ctg\t",$i+1,"\t.\t",substr($seq,$i,1),"\t",substr($seq,$i,1),"$sub\t*\t*\t*\t*\t1:1:1:0:0:10:10:0\n";
      }
    $i+=$size+1;
    }
  }
}
}
          
          

