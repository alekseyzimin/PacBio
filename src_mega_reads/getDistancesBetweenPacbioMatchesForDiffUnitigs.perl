#!/usr/bin/perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/perl
for ($i=0; $i<=$#ARGV; ++$i) {
    if ($ARGV[$i] eq "-l") {
	++$i;
	$lengthFile = $ARGV[$i];
	next; }
    elsif ($ARGV[$i] eq "-p") {
	++$i;
	$posmapFile = $ARGV[$i];
	next; }
    elsif ($ARGV[$i] eq "-test") {
	$workDir = "/genome3/raid/alekseyz/PB_ScerW303/mega-reads/assembly";
	$lengthFile = "$workDir/assembly.CA-bogart-bogart.200.1.unitig_lengths.txt";
	$posmapFile = "$workDir/assembly.CA-bogart-bogart.200.1.posmap.txt"; 
	next; }
    print STDERR "The args are incorrect\n";
    exit (1);
}

open (FILE, $lengthFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $length[$flds[0]] = $flds[1];
}
close (FILE);

open (FILE, $posmapFile);
while ($line = <FILE>) {
    next unless ($line =~ /^m/);
    chomp ($line);
    @flds = split (" ", $line);
    ($pacbio, $part) = ($flds[0] =~ /^([^\.]+)\.(.+)$/);
    $unitig = $flds[1];
    $begin = $flds[2];
    $end = $flds[3];
    $lines{$pacbio} .= "$line\n";
    if (! $unitig{$pacbio}) {
	($pacbioOffset1, $pacbioOffset2) = ($pacbio =~ /\D(\d+)_(\d+)$/);
	$pacbioLength{$pacbio} = $pacbioOffset2 - $pacbioOffset1;
	$unitig{$pacbio} = $unitig; }
    elsif ($unitig{$pacbio} != $unitig) {
	$isNeeded{$pacbio} = 1; }
}
close (FILE);

@pacbios = keys %isNeeded;

for (@pacbios) {
    $pacbio = $_;
    $tline = $lines{$pacbio};
    $pacbioLength = $pacbioLength{$pacbio};
#    print "tline = \n";
#    print "$tline\n"; # For printing out as before
#    next;
    chomp ($tline);
    @lines = split (/\n/, $tline);
    @lines2 = @lines3 = ();
    for (@lines) {
	$line = $_;
	($piece) = ($line =~ /^\S+\.(\d+)_\d+\s/);
	$lines3[$piece] = "$line\n"; }
    for ($i=0; $i<=$#lines3; ++$i) {
	if ($lines3[$i] =~ /\S/) {
	    push (@lines2, $lines3[$i]); }
    }
    # For first match to unitig
    $line1 = $lines2[0];
    chomp ($line1);
    @flds1 = split (" ", $line1);
    if ($flds1[2] < $flds1[3]) { # Ori is 'F'
	$sep1 = $length[$flds1[1]] - $flds1[2]; }
    else {
	$sep1 = $flds1[2]; }
    # For last match to unitig
    $line2 = $lines2[-1];
    chomp ($line2);
    @flds2 = split (" ", $line2);
    $lastUnitig = $flds2[1];
    if ($flds2[2] < $flds2[3]) { # Ori is 'F'
	$sep2 = $flds2[3]; }
    else {
	$sep2 = $length[$flds2[1]] - $flds2[3]; }
    $lastUnitig = $flds2[1];
    # Only add unitig lengths until we see the last consecutive matches to
    # the last unitig
    for ($lastIndexToConsider=$#lines2; $lastIndexToConsider>0; --$lastIndexToConsider) {
	$line = $lines2[$lastIndexToConsider-1];
	chomp ($line);
	@flds = split (" ", $line);
	last if ($flds[1] != $lastUnitig);
    }
    $midUniLengthSum = 0;
    for ($i=0; $i<$lastIndexToConsider-1; ++$i) {
	$line1 = $lines2[$i];
	$line2 = $lines2[$i+1];
	chomp ($line1);
	chomp ($line2);
	@flds1 = split (" ", $line1);
	@flds2 = split (" ", $line2);
	next if ($flds1[1] == $flds2[1]); # Same unitig
	$midUniLengthSum += $length[$flds2[1]];
#	print "uni $flds2[1] $length[$flds2[1]] $midUniLengthSum\n";
    }
    $sumImpliedLengths = $midUniLengthSum + $sep1 + $sep2;
    $impliedLengthToPacbioLength = $sumImpliedLengths / $pacbioLength;
    
    print "$pacbio $sumImpliedLengths $pacbioLength $impliedLengthToPacbioLength\n";
#    print @lines2, "\n";
}

