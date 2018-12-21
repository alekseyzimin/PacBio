#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
#
$last_pb_read="";
$last_pb_start=-1;
$last_seq="";
$max_read_length=400;
@mate_pairs=();
while($line=<STDIN>){
    chomp($line);
    if($line=~/^\>/){
($readname)=split(/\s+/,substr($line,1));
$line=<STDIN>;
chomp($line);
my $seq=$line;
my ($pb_read,$pb_coords)=split(/\./,$readname);
my  ($pb_start,$pb_len)=split(/\_/,$pb_coords);
if($pb_read eq $last_pb_read){#we have a mate pair
    my $len1=$max_read_length;
    my $len2=$max_read_length;
    $len1=int(length($last_seq)*.9) if(length($last_seq)<$len1);
    $len2=int(length($seq)*.9) if(length($seq)<$len2);
    my $gap_len=$pb_start-($last_pb_start+$last_pb_len);
    my @ttt=($pb_read,substr($last_seq,20,$len1),substr(reverse_complement($seq),20,$len2),$pb_len+$last_pb_len+$gap_len*.5);
    push(@mate_pairs,\@ttt);
}
$last_pb_read=$pb_read;
$last_pb_start=$pb_start;
$last_pb_len=$pb_len;
$last_seq=$seq;
    }
}

my $max_length=0;
foreach $mp(@mate_pairs){
	$max_length=$$mp[3] if($$mp[3]>$max_length);
}
$max_length+=500;

print STDOUT "{VER\n";
print STDOUT "ver:2\n";
print STDOUT "}\n";

for($i=2000;$i<=$max_length;$i+=1000){
    print STDOUT "{LIB\n";
    print STDOUT "act:A\n";
    print STDOUT "acc:mr_",int(($i+500)/1000),"\n";
print STDOUT "ori:I\n";
print STDOUT "mea:$i\n";
print STDOUT "std:",int($i*.075),"\n";
print STDOUT "src:\n";
print STDOUT ".\n";
print STDOUT "nft:4\n";
print STDOUT "fea:\n";
print STDOUT "doRemoveChimericReads=1\n";
print STDOUT "doRemoveSpurReads=1\n";
print STDOUT "isNotRandom=1\n";
print STDOUT "constantInsertSize=1\n";
print STDOUT ".\n";
print STDOUT "}\n";
}

$mp_num=-1;
foreach $mp(@mate_pairs){
    next if($$mp[3]<2500);
$mp_num++;
my $readname1=$$mp[0].".".$mp_num."F";
my $readname2=$$mp[0].".".$mp_num."R";
my $libId="mr_".int(($$mp[3]+500)/1000);
my $sequence1=$$mp[1];
my $sequence2=$$mp[2];
next if(substr($sequence1,0,250) eq substr($sequence2,0,250));

#print join(" ",@$mp),"\n";
print STDOUT "{FRG\n";
print STDOUT "act:A\n";
print STDOUT "acc:$readname1\n";
print STDOUT "rnd:0\n";
print STDOUT "sta:G\n";
print STDOUT "lib:$libId\n";
print STDOUT "pla:0\n";
print STDOUT "loc:0\n";
print STDOUT "src:\n.\n";
print STDOUT "seq:\n$sequence1\n.\n";
$sequence1 =~ tr/ACGTNacgtn/XXXXXDDDDD/;# create fake quality scores
print STDOUT "qlt:\n$sequence1\n.\n";
print STDOUT "hps:\n.\n";
print STDOUT "clv:0,",length($sequence1),"\n";
print STDOUT "clr:0,",length($sequence1),"\n";
print STDOUT "}\n";
print STDOUT "{FRG\n";
print STDOUT "act:A\n";
print STDOUT "acc:$readname2\n";
print STDOUT "rnd:0\n";
print STDOUT "sta:G\n";
print STDOUT "lib:$libId\n";
print STDOUT "pla:0\n";
print STDOUT "loc:0\n";
print STDOUT "src:\n.\n";
print STDOUT "seq:\n$sequence2\n.\n";
$sequence2 =~ tr/ACGTNacgtn/XXXXXDDDDD/;# create fake quality scores
print STDOUT "qlt:\n$sequence2\n.\n";
print STDOUT "hps:\n.\n";
print STDOUT "clv:0,",length($sequence2),"\n";
print STDOUT "clr:0,",length($sequence2),"\n";
print STDOUT "}\n";
print STDOUT "{LKG\n";
print STDOUT "act:A\n";
print STDOUT "frg:$readname1\n";
print STDOUT "frg:$readname2\n";
print STDOUT "}\n";
}

sub reverse_complement{
    my $sequence=$_[0];
    $sequence=reverse($sequence);
    $sequence=~tr/ACGTNacgtn/TGCANtgcan/;
    return($sequence);
}

