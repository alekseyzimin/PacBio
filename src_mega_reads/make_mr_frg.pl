#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl

my $libId=$ARGV[0];
my $min_len_output=64;
$min_len_output=$ARGV[1] if($#ARGV>0);

print STDOUT "{VER\n";
print STDOUT "ver:2\n";
print STDOUT "}\n";
print STDOUT "{LIB\n";
print STDOUT "act:A\n";
print STDOUT "acc:$libId\n";
print STDOUT "ori:I\n";
print STDOUT "mea:3000\n";
print STDOUT "std:300\n";
print STDOUT "src:\n";
print STDOUT ".\n";
print STDOUT "nft:3\n";
print STDOUT "fea:\n";
print STDOUT "doTrim_initialNone=1\n";
print STDOUT "doRemoveChimericReads=1\n";
print STDOUT "doRemoveSpurReads=1\n";
print STDOUT ".\n";
print STDOUT "}\n";

while($line1=<STDIN>)
{
chomp($line1);
if($line1 =~ /^>/)
{
$header=substr($line1,1);
@f=split(/\s+/,$header);
$readname1=$f[0];
$line1=<STDIN>;
chomp($line1);
$sequence1=$line1;
$clr1=0;
$clr2=length($sequence1);
next if(length($sequence1)<$min_len_output);

        print STDOUT "{FRG\n";
        print STDOUT "act:A\n";
        print STDOUT "acc:$readname1\n";
        print STDOUT "rnd:1\n";
        print STDOUT "sta:G\n";
        print STDOUT "lib:$libId\n";
        print STDOUT "pla:0\n";
        print STDOUT "loc:0\n";
        print STDOUT "src:\n.\n";
        print STDOUT "seq:\n$sequence1\n.\n";
$sequence1 =~ tr/ACGTNacgtn/XXXXXDDDDD/;# create fake quality scores
        print STDOUT "qlt:\n$sequence1\n.\n";
        print STDOUT "hps:\n.\n";
        print STDOUT "clv:$clr1,$clr2\n";
        print STDOUT "clr:$clr1,$clr2\n";
        print STDOUT "}\n";

}
}
