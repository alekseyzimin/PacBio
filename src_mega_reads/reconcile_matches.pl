#!/usr/bin/env perl


my $tol_factor=1.2;#tolerance factor around gap size
my $tol_min=200;#minimum tolerance for small gap
my $scf="";
my $line=<>,@l=();
my @f=split(/\s+/,$line);
my $scf=$f[-2];
push(@l,$line);

while($line=<>){
$line=~s/^\s+//;
@f=split(/\s+/,$line);
if(not($f[-2] eq $scf)){
process_lines(@l);
@l=();
$scf=$f[-2];
}
push(@l,$line);
}
process_lines(@l);

sub process_lines{
my @lines=@_;
my $gap_before=10000;
my $gap_after=3000;
my @l1,@l2,@l3;
@l2=split(/\s+/,$lines[0]);
@l3=split(/\s+/,$lines[1]);
$gap_after=($l3[0]-$l2[1])*$tol_factor;
$gap_after=$tol_min if($gap_after<$tol_min);
output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
for($i=1;$i<$#lines;$i++){
@l1=split(/\s+/,$lines[$i-1]);
@l2=split(/\s+/,$lines[$i]);
@l3=split(/\s+/,$lines[$i+1]);

$gap_before=($l2[0]-$l1[1])*$tol_factor;
$gap_before=$tol_min if($gap_before<$tol_min);

$gap_after=($l3[0]-$l2[1])*$tol_factor;
$gap_after=$tol_min if($gap_after<$tol_min);

output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],1);
}

@l1=split(/\s+/,$lines[$#lines-1]);
@l2=split(/\s+/,$lines[$#lines]);

$gap_before=($l2[0]-$l1[1])*$tol_factor;
$gap_before=$tol_min if($gap_before<$tol_min);
$gap_after=10000;
output_coords($gap_before,$gap_after,$l2[3],$l2[4],$l2[12],$l2[-2],$l2[-1],0);
}

sub output_coords{
my $start,$end,$dir;
my ($gap_b,$gap_a,$s,$e,$len,$scf,$ctg)=@_;
if($s<$e){#forward match
$dir="f";
if($s<$gap_b){
$start=0;
$gap_b-=($s-1);
}else{
$start=$s;
}
if(($len-$e)<$gap_a){
$end=$len;
$gap_a-=($len-$e);
}else{
$end=$e;
}
}else{
$dir="r";
if($e<$gap_a){
$start=0;
$gap-=($e-1);
}else{
$start=$e;
}
if(($len-$s)<$gap_b){
$end=$len;
$gap_b-=($len-$s);
}else{
$end=$s;
}
}
print "$scf $ctg $start $end $dir $gap_b $gap_a $len\n";
}
