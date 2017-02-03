#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
# $file = "/genome3/raid/alekseyz/PB_ScerW303/mega-reads-new/mrN6_max_2p_1.70.13.25.0.1.allowed";
$errorRateAllowed = .1;
$errorMin = 10;
$correctnessCodeForSingletons = 1;
&processArgs;

#open (FILE, $file);
&initializeData;
$line = <STDIN>;
&loadDataFromLine;
$indexHold = $index;
&addDataToRecords;
while ($line = <STDIN>) {
    &loadDataFromLine;
    if ($index ne $indexHold) {
	&analyzeGroup;
	&initializeData;
	$indexHold = $index; }
    &addDataToRecords;
}
&analyzeGroup;

sub findStartingNumber # Calculates $gapMin, $gapMax
{
    my ($mostInUnion, $i, $j, $currentInUnion, $absGap, $allowedGapDiff);
    my ($currentGapMin, $currentGapMax, $gapMax, $gapMin, $bestGapSize);
    my ($gapSum, $mean);

    $mostInUnion = 0;
    for ($i=0; $i<=$#gaps; ++$i) {
	$currentInUnion = 0;
	$absGap = abs ($gaps[$i]);
	$allowedGapDiff = $errorMin;
	if ($errorRateAllowed * $absGap > $errorMin) {
	    $allowedGapDiff = $errorRateAllowed * $absGap; }
	$allowedGapDiff*=2;  #because there may be a center that satisfies the condition for two numbers
	for ($j=$i; $j>=0; --$j) {
	    last if ($gaps[$j] < $gaps[$i] - $allowedGapDiff);
	    $currentGapMin = $gaps[$j];
	    $currentGapMax = $gaps[$j];
	    ++$currentInUnion; }
	for ($j=$i+1; $j<=$#gaps; ++$j) {
	    last if ($gaps[$j] > $gaps[$i] + $allowedGapDiff);
	    $currentGapMax = $gaps[$j];
	    ++$currentInUnion; }
	if ($currentInUnion > $mostInUnion) {
	    $bestGapSize = $gaps[$i];
	    $mostInUnion = $currentInUnion;
	    $gapMin = $currentGapMin;
	    $gapMax = $currentGapMax; }
    }

    # Now we calculate the average of the bin and recalculate the bin using
    # the mean as the center
    for ($i=0; $i<=$#gaps; ++$i) {
	next if ($gaps[$i] < $gapMin);
	last if ($gaps[$i] > $gapMax);
	$gapSum += $gaps[$i]; }
    $mean = $gapSum / $mostInUnion;
#     print "gapMin = $gapMin gapMax = $gapMax mean = $mean\n";
    $absGap = abs ($mean);
    $allowedGapDiff = $errorMin;
    if ($errorRateAllowed * $absGap > $errorMin) {
	$allowedGapDiff = $errorRateAllowed * $absGap; }
    $gapMin = -1000000;
    $mostInUnion = 0;
    for ($i=0; $i<=$#gaps; ++$i) {
	next if ($gaps[$i] < $mean - $allowedGapDiff);
	last if ($gaps[$i] > $mean + $allowedGapDiff);
	if ($gapMin <= -1000000) {
	    $gapMin = $gaps[$i]; }
	$gapMax = $gaps[$i];
	++$mostInUnion;
    }
    # End of section recomputing the average of the bin and recalculating

    $minGood = $gapMin;
    $maxGood = $gapMax;
    if ($mostInUnion == 1) {
	return ("bad"); }
    return ("good");
}

sub analyzeGroup
{
    my ($line, @flds, $goodCode);
    @gaps = sort byNum @gaps;
    if ($#gaps == 0) {
	@flds = split (" ", $lines[0]);
	print "@flds[0..3] $correctnessCodeForSingletons\n";
	return; }
    $goodCode = &findStartingNumber; # Sets $minGood, $maxGood
    if ($goodCode eq "bad") {
	for (@lines) {
	    $line = $_;
	    @flds = split (" ", $line);
	    print "@flds[0..3] 0\n"; }
	return;
    }
    for (@lines) {
	$line = $_;
	@flds = split (" ", $line);
	print "@flds[0..3] ";
	if (($flds[1] >= $minGood) && ($flds[1] <= $maxGood)) {
	    print "1\n"; }
	else {
	    print "0\n"; }
    }
}

sub initializeData
{
    @gaps = @lines = ();
}

sub loadDataFromLine
{
    my (@flds);
    chomp ($line);
    @flds = split (" ", $line);
    $index = "@flds[2..3]";
    $sep = $flds[1];
    $accepted = $flds[6];
}

sub addDataToRecords
{
    push (@gaps, $sep);
    push (@lines, $line);
}

sub byNum
{
    return ($a <=> $b);
}

sub processArgs
{
    for ($i=0; $i<=$#ARGV; ++$i) {
	if ($ARGV[$i] =~ /^[-]/) {
	    if ($ARGV[$i] =~ /^[-][-]?[hH](elp)?$/) {
		&reportUsage; }
	    if ($i == $#ARGV) {
		print STDERR "You cannot end the args with a flag. Bye!\n";
		&reportUsage; }
#	    if (($ARGV[$i] eq "-f") || ($ARGV[$i] eq "--file")) {
#		++$i;
#		$file = $ARGV[$i];
#		next; }
	    if ($ARGV[$i] eq "--min-range-radius") {
		++$i;
		$errorMin = $ARGV[$i];
		next; }
	    if ($ARGV[$i] eq  "--min-range-proportion") {
		++$i;
		$errorRateAllowed = $ARGV[$i];
		next; }
	    else {
		print STDERR "Unknown flag: $ARGV[$i]. Bye!\n";
		&reportUsage; }
	}
    }
#    if ($file !~ /\S/) {
#	print STDERR "File not specified. Bye!\n";
#	&reportUsage; }
#    if (! -e $file) {
#	print STDERR "File \"$file\" doesn't exist. Bye!\n";
#	&reportUsage; }
#    if (-s $file == 0) {
#	print STDERR "File \"$file\" has size 0. Bye!\n";
#	&reportUsage; }
    if ($errorMin !~ /^\d+$/) {
	print STDERR "The flag --min-range-radius must have a non-negative integer argument. Given value is ${errorMin}. Bye!\n";
	exit (1); }
    if ($errorRateAllowed !~ /^0*\.[0-9]*[1-9][0-9]*$/) {
	print STDERR "The flag --min-range-proportion must have a decimal argument that is less than 1 and greater than 0. Bye!\n";
	exit (1); }
}

sub reportUsage
{
    print STDERR "To run this file you must specify the following arg:\n  -f inputFile   (inputFile is the name of the file with k-unitig pairs and gap sizes)\n  --min-range-radius int  (int >= 0, specifies dist from center that gap sizes are accepted) (default: 10)\n  --min-range-proportion float  (0<float<1, specifies proportion of center gap size that we allow for accepted gap sizes) (default: .1)\n  -h or --help is help\n";
    exit (1);
}
