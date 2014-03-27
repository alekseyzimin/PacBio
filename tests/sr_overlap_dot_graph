#!/usr/bin/env perl

use strict;
use warnings;

my (@minOffsets, @maxOffsets, @numKMers, @superReadLengths, @beginSupers, @endSupers, @superReadNames, @mersCounts);
my $pacbioName;

while(my $line = <STDIN>) {
    chomp ($line);
    my @flds = split (" ", $line);
    push (@minOffsets, $flds[0]);
    push (@maxOffsets, $flds[1]);
    push (@numKMers, $flds[4]);
    push (@superReadLengths, $flds[10]);
    if (defined ($pacbioName)) {
	if ($flds[11] ne $pacbioName) {
	    print "You screwed up! Only use for ONE (1) pacbio read at a time. Bye!\n";
	    exit (1); }
    }
    else {
	$pacbioName = $flds[11]; }
    my $superReadName = $flds[12];
    if ($flds[2] > $flds[3]) {
	$superReadName = &reverseSuperReadName ($superReadName); 
	$flds[2] = $flds[10] + 1 - $flds[2];
	$flds[3] = $flds[10] + 1 - $flds[3];
    }
    push (@beginSupers, $flds[2]);
    push (@endSupers, $flds[3]);
    push (@superReadNames, $superReadName);
    splice(@flds, 0, 13);
    push (@mersCounts, \@flds);
}

print "digraph \"$pacbioName\" {\n";
print "node [fontsize=10];\n";
for (my $i=0; $i< @superReadNames; ++$i) {
#    print "$i $minOffsets[$i] $maxOffsets[$i] $superReadNames[$i]\n";
    print "$i [tooltip=\"$superReadNames[$i]\",label=\"($i) L$superReadLengths[$i] #$numKMers[$i]\\nP ($minOffsets[$i],$maxOffsets[$i])\\nS ($beginSupers[$i],$endSupers[$i])\"];\n";
}

for (my $i=0; $i< @superReadNames; ++$i) {
    my $minOffset = $minOffsets[$i];
    my $maxOffset = $maxOffsets[$i];
    my $super = $superReadNames[$i];
    my @tflds = split (/_/, $super);
    my $lastKUni = $tflds[-1];
    for (my $j=$i+1; $j< @superReadNames; ++$j) {
	last if ($minOffsets[$j] >= $maxOffset);
	my $localSuperReadName = $superReadNames[$j];
	my $index = index ($localSuperReadName, $lastKUni);
	next unless ($index >= 0);
	my $lenOfSubstring = $index + length ($lastKUni);
	my $lastToFind = substr ($super, -$lenOfSubstring);
	$index = index ($localSuperReadName, $lastToFind);
	next unless ($index == 0);
        my @unitigs = split(/_/, $lastToFind); # nb k-unitigs in common between nodes
        my $mers_shared = 0;
        my $mers_info = $mersCounts[$j];
        for(my $k = 0; $k < @unitigs; ++$k) {
          $mers_shared += $$mers_info[2 * $k];
          $mers_shared -= $$mers_info[2 * $k - 1] if $k > 0;
        }
	print "$i -> ${j} [tooltip=\"$lastToFind\", label=\"$mers_shared\"];\n";
    }
}
print "}\n";
	

sub reverseSuperReadName
{
    my ($super) = @_;

    my @flds = split (/_/, $super);
    my $superOut = "";
    for (my $i=$#flds; $i>=0; --$i) {
        my $fld = $flds[$i];
	$fld =~ tr/FR/RF/;
	if ($i != $#flds) {
	    $superOut .= "_"; }
	$superOut .= $fld; }
    return ($superOut);
}

