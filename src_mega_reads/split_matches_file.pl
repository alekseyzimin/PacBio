#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
#This file splits the pacbio matches file into chunks for parallel processing
my $chunksize=$ARGV[0];
my $prefix=$ARGV[1];

my $counter=0;
my $file_counter=0;
open(FILE,">$prefix.$file_counter");
while($line=<STDIN>){
    if($line =~ /^>/){
	if($counter>$chunksize){
	$file_counter++;
	open(FILE,">$prefix.$file_counter");
	$counter=0;
	}else{
	$counter++;
	}
}
print FILE $line;
}

