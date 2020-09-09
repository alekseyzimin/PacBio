#!/usr/bin/env perl
#this code converts the output of the split scaffolds to reconciled.txt file format
$current_chr="";
@lines=();
while($line1=<STDIN>){
chomp($line1);
($chr,$rest)=split(/:/,$line1);
if(not($current_chr eq $chr)){
  $current_chr=$chr;
  if(scalar(@lines)>0){
    $line=$lines[0];
    ($ctg,$size)=split(/\s+/,$line);
    ($chr,$coords)=split(/:/,$ctg);
    @coords=split(/-/,$coords);
    $gap_before=100;
    if(scalar(@lines)>1){
        $linea=$lines[1];
        ($ctga,$sizea)=split(/\s+/,$linea);
        ($chra,$coordsa)=split(/:/,$ctga);
        @coordsa=split(/-/,$coordsa);
        $gap_after=$coordsa[0]-$coords[1]-1;
        print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
        for($i=1;$i<=$#lines-1;$i++){
          $gap_before=$gap_after;
          $line=$lines[$i];
          ($ctg,$size)=split(/\s+/,$line);
          ($chr,$coords)=split(/:/,$ctg);
          @coords=split(/-/,$coords);
          $linea=$lines[$i+1];
          ($ctga,$sizea)=split(/\s+/,$linea);
          ($chra,$coordsa)=split(/:/,$ctga);
          @coordsa=split(/-/,$coordsa);
          $gap_after=$coordsa[0]-$coords[1]-1;
          print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
        }
        $gap_before=$gap_after;
        $gap_after=100;
        $line=$lines[$#lines];
        ($ctg,$size)=split(/\s+/,$line);
        print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
    }else{
        $gap_after=100;
        print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
    }
    }
    @lines=();
  } 
  push(@lines,$line1);
}
#last one
  if(scalar(@lines)>0){
    $line=$lines[0];
    ($ctg,$size)=split(/\s+/,$line);
    ($chr,$coords)=split(/:/,$ctg);
    @coords=split(/-/,$coords);
    $gap_before=100;
    if(scalar(@lines)>1){
        $linea=$lines[1];
        ($ctga,$sizea)=split(/\s+/,$linea);
        ($chra,$coordsa)=split(/:/,$ctga);
        @coordsa=split(/-/,$coordsa);
        $gap_after=$coordsa[0]-$coords[1]-1;
        print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
        for($i=1;$i<=$#lines-1;$i++){
          $gap_before=$gap_after;
          $line=$lines[$i];
          ($ctg,$size)=split(/\s+/,$line);
          ($chr,$coords)=split(/:/,$ctg);
          @coords=split(/-/,$coords);
          $linea=$lines[$i+1];
          ($ctga,$sizea)=split(/\s+/,$linea);
          ($chra,$coordsa)=split(/:/,$ctga);
          @coordsa=split(/-/,$coordsa);
          $gap_after=$coordsa[0]-$coords[1]-1;
          print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
        }
        $gap_before=$gap_after;
        $gap_after=100;
        $line=$lines[$#lines];
        ($ctg,$size)=split(/\s+/,$line);
        print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
    }else{
        $gap_after=100;
        print "$chr $ctg 1 $size f $gap_before $gap_after $size\n";
    }
    }

