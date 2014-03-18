#!/usr/bin/env perl
while ($line = <STDIN>) {
    chomp ($line);
    @flds = split (" ", $line);
    push (@minOffsets, $flds[0]);
    push (@maxOffsets, $flds[1]);
    push (@numKMers, $flds[4]);
    push (@superReadLengths, $flds[10]);
    if (defined ($pacbioName)) {
	if ($flds[-2] ne $pacbioName) {
	    print "You screwed up! Only use for ONE (1) pacbio read at a time. Bye!\n";
	    exit (1); }
    }
    else {
	$pacbioName = $flds[-2]; }
    $superReadName = $flds[-1];
    if ($flds[2] > $flds[3]) {
	$superReadName = &reverseSuperReadName ($superReadName); 
	$flds[2] = $flds[10] + 1 - $flds[2];
	$flds[3] = $flds[10] + 1 - $flds[3];
    }
    push (@beginSupers, $flds[2]);
    push (@endSupers, $flds[3]);
    push (@superReadNames, $superReadName);
}

print "digraph \"$pacbioName\" {\n";
print "node [fontsize=10];\n";
for ($i=0; $i<=$#superReadNames; ++$i) {
#    print "$i $minOffsets[$i] $maxOffsets[$i] $superReadNames[$i]\n";
    print "$i [tooltip=\"$superReadNames[$i]\",label=\"($i) L$superReadLengths[$i] #$numKMers[$i]\\nP ($minOffsets[$i],$maxOffsets[$i])\\nS ($beginSupers[$i],$endSupers[$i])\"];\n";
}

for ($i=0; $i<=$#superReadNames; ++$i) {
    $minOffset = $minOffsets[$i];
    $maxOffset = $maxOffsets[$i];
    $super = $superReadNames[$i];
    @tflds = split (/_/, $super);
    $lastKUni = $tflds[-1];
    for ($j=$i+1; $j<=$#superReadNames; ++$j) {
	last if ($minOffsets[$j] >= $maxOffset);
	$localSuperReadName = $superReadNames[$j];
	$index = index ($localSuperReadName, $lastKUni);
	next unless ($index >= 0);
	$lenOfSubstring = $index + length ($lastKUni);
	$lastToFind = substr ($super, -$lenOfSubstring);
	$index = index ($localSuperReadName, $lastToFind);
	next unless ($index == 0);
	print "$i -> ${j} [tooltip=\"$lastToFind\"];\n";
    }
}
print "}\n";
	

sub reverseSuperReadName
{
    my ($super) = @_;
    my (@flds, $superOut, $i, $fld);

    @flds = split (/_/, $super);
    $superOut = "";
    for ($i=$#flds; $i>=0; --$i) {
	$fld = $flds[$i];
	$fld =~ tr/FR/RF/;
	if ($i != $#flds) {
	    $superOut .= "_"; }
	$superOut .= $fld; }
    return ($superOut);
}

