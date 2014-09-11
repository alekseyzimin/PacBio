#!/usr/bin/perl
$fn = "/genome3/raid/alekseyz/PB_ScerW303/mega-reads/coords_m300_k15_10x_70_B15_U1_mm";
$minOvlForMatchingRgns = 100;
$minOvl = 0;
$minMatchLenForImpliedMatch = 25;
open (FILE, $fn);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    if ($line =~ /^>/) {
	if ($pacbio) {
	    if ($isFirstLine) {
#		print "$pacbio has no matches to super-reads\n";
	    }
	    else {
		# Must report if nothing at end
		$diff = $end - $pacbioLen;
		if ($diff < 0) {
#		    print "$pacbio $diff 0\n";
		}
		&reportNonOverlappedGaps;
	    }
	}
	($pacbio) = ($line =~ /^.\d+\s+(\S.+)$/);
	$begin = $end = 0;
	@beginsToCover = @endsToCover = ();
	@impliedBegins = @impliedEnds = ();
	$isFirstLine = 1;
	next; }
    if ($isFirstLine) {
	$pacbioLen = $flds[9];
	$isFirstLine = 0; }

    $begin = $flds[0]-1;
    if ($flds[1]-$flds[0]+1 >= $minMatchLenForImpliedMatch) {
	$impliedBegin = $flds[0] - $flds[2];
	if ($impliedBegin < 0) { $impliedBegin = 0; }
	$impliedEnd = $flds[1] + ($flds[10] - $flds[3]);
	if ($impliedEnd > $pacbioLen) { $impliedEnd = $pacbioLen; }
	push (@impliedBegins, $impliedBegin);
	push (@impliedEnds, $impliedEnd);
#       print "impliedBegin = $impliedBegin impliedEnd = $impliedEnd\n";
    }
    if ($begin>$end-$minOvlForMatchingRgns) {
	if ($end > 0) {
	    if ($end < $begin) {
		push (@beginsToCover, $end);
		push (@endsToCover, $begin); }
	    else {
		push (@beginsToCover, $begin);
		push (@endsToCover, $end); }
	}
#	print "$pacbio $end $begin\n";
    }
    if ($end < $flds[1]) {
	$end = $flds[1]; }
}

if ($pacbio) {
    if ($isFirstLine) {
#	print "$pacbio has no matches to super-reads\n";
    }
    else {
	# Must report if nothing at end
	$diff = $end - $pacbioLen;
	if ($diff < 0) {
#	    print "$pacbio $diff 0\n";
	}
	&reportNonOverlappedGaps;
    }
}

sub reportNonOverlappedGaps
{
    my ($i, $j, @indices, @reducedImpliedBegins, @reducedImpliedEnds, $intervalBegin, $intervalEnd);

    return if ($#impliedBegins < 0);
    return if ($#beginsToCover < 0);
    @indices = (0..$#impliedBegins);
    @indices = sort spcl @indices;
    @reducedImpliedBegins = @reducedImpliedEnds = ();
    @reducedImpliedBegins = ($impliedBegins[$indices[0]]);
    @reducedImpliedEnds = ($impliedEnds[$indices[0]]);
    for ($i=1; $i<=$#impliedBegins; ++$i) {
	next if ($impliedEnds[$indices[$i]] <= $reducedImpliedEnds[-1]);
	push (@reducedImpliedBegins, $impliedBegins[$indices[$i]]);
	push (@reducedImpliedEnds, $impliedEnds[$indices[$i]]);
    }
    $i = $j = 0;
    $intervalBegin = -1;
    $intervalEnd = -1;
    while (1) {
	if ($reducedImpliedEnds[$j] < $endsToCover[$i] + $minOvl) {
	    ++$j;
	    last if ($j > $#reducedImpliedEnds);
	    next; }
	if ($reducedImpliedBegins[$j] <= $beginsToCover[$i] - $minOvl) {
	    if ($intervalBegin < 0) {
		$intervalBegin = $beginsToCover[$i];
		$intervalEnd = $endsToCover[$i]; }
#	    print "$pacbio $beginsToCover[$i] $endsToCover[$i] $reducedImpliedBegins[$j] $reducedImpliedEnds[$j]\n";
#	    print "intervalBegin = $intervalBegin\n";
	    if ($beginsToCover[$i] > $intervalEnd) {
		print "BREAK $pacbio $intervalBegin $intervalEnd\n";
		$intervalBegin = $beginsToCover[$i]; }
	    if ($endsToCover[$i] > $intervalEnd) {
		$intervalEnd = $endsToCover[$i]; }
	}
	++$i;
	last if ($i > $#beginsToCover);
    }
    if ($intervalBegin >= 0) {
	print "BREAK $pacbio $intervalBegin $intervalEnd\n"; }
}

sub spcl
{
    if ($impliedBegins[$a] != $impliedBegins[$b]) {
	return ($impliedBegins[$a] <=> $impliedBegins[$b]); }
    return ($impliedEnds[$b] <=> $impliedEnds[$a]);
}

