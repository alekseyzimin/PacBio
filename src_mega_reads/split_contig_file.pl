#!/usr/bin/perl
# Pass the name of a fasta file as an arg and this gives the number of
# bases in each read of the fasta file
open (FILE, "$ARGV[1]");
@f=split(/\//,$ARGV[1]);
my $outfilename="$ARGV[0]/$f[-1]";
my $fileindex=1;
my $basestodate=0;
my $basestosplit=$ARGV[2];
my $currentfilename="$outfilename\.$fileindex";
open(OUTFILE,">$currentfilename");
while ($line = <FILE>) 
{
    	if ($line =~ /^>/) 
	{	
		#check if need to open new file
		if($basestodate>$basestosplit)
		{
		print "$basestodate bases output in file ",$currentfilename,"\n";
		$fileindex++;
		close(OUTFILE);
		$currentfilename="$outfilename\.$fileindex";
		open(OUTFILE,">$currentfilename");
		$basestodate=0;
		}
    	}
    	else 
	{
		$len = length ($line)-1;
		$basestodate +=$len;
    	}
	print OUTFILE $line;
}

print "$basestodate bases output in file ",$currentfilename,"\n";
close(OUTFILE);
