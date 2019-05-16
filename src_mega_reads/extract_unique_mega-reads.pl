#!/usr/bin/env perl
#this code extracts the unique mega-reads from the output of create_mega_reads, normalizes their names and prints them out on STDOUT, lengths are printed out on STDERR
while($line=<STDIN>){
  chomp($line);
  if(substr($line,0,1) eq ">"){
    $pb=substr($line,1);
  }else{
    ($junk0,$junk1,$junk2,$junk3,$junk4,$junk5,$junk6,$junk7,$mega_read,$junk9,$sequence)=split(/\s+/,$line);
    @kunis=split(/_/,$mega_read);
    if(substr($kunis[0],0,-1)>substr($kunis[-1],0,-1)){
      $mega_read=join("_",reverse(@kunis));
      $mega_read=~tr/FR/RF/;
      $sequence=reverse($sequence);
      $sequence=~tr/ACGTNacgtn/TGCANtgcan/;
    }
    if(length($mega_read)<length($sequence)){
      $mega_read_index=$mega_read;
    }else{
      $mega_read_index=$sequence;
    }
    if(not(defined($out{$mega_read_index}))){
      print STDOUT ">$mega_read\n$sequence\n";
      print STDERR "$mega_read ",length($sequence),"\n";
      $out{$mega_read_index}=1;
    }
  }
}
