#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
$errorRateAllowed = .1;
$errorMin = 10;
$correctnessCodeForSingletons = 0;
&processArgs;
my %groups;

#load data
while($line=<STDIN>){
  chomp($line);
  @f=split(/\s+/,$line);
  push(@{$groups{"$f[2] $f[3]"}},$line);
}

foreach $group(keys %groups){
  my $groupSize=scalar(@{$groups{$group}});
  if($groupSize==1){
    print ${$groups{$group}}[0]," $correctnessCodeForSingletons\n";
  }elsif($groupSize==2){
    my $group_center=(${$groups{$group}}[0]+${$groups{$group}}[1])/2;
    $group_center=0.00001 if($group_center==0);
    my $groupCode=0;
    $groupCode=1 if(abs(${$groups{$group}}[0]-$group_center)<=$errorMin*2 || abs((${$groups{$group}}[0]-$group_center)/$group_center)<=$errorRateAllowed*2);
    foreach $l(@{$groups{$group}}){
      print $l," $groupCode\n";
    }
  }else{#more than three elements:  sort, find the median then look for radius of at least $errorMin or median*$errorRateAllowed around the median
    @lines_sorted=sort byGap @{$groups{$group}};
    @line_gaps=();
    foreach $l(@lines_sorted){
      my @f=split(/\s+/,$l);
      push(@line_gaps,$f[1]);
    }

    my $new_median_value=$line_gaps[int(scalar(@lines_sorted)/2)];
    $new_median_value+=0.000001 if($new_median_value==0);
    my $median_value=100000;
    my $exit_code=0;
    my $radius=0;
    my $max_iterations=5;
    my $iteration=0;
    #print "$new_median_value $median_value $iteration\n";
    while(abs(($median_value-$new_median_value)/$new_median_value)>$errorRateAllowed && abs($median_value-$new_median_value)>$errorMin && $iteration<$max_iterations ){
      $iteration++;
      #print "DEBUG $iteration $new_median_value\n";
      my @line_gaps_new=();
      $median_value=$new_median_value;
      $radius=abs($median_value*$errorRateAllowed);
      $radius=$errorMin if($radius<$errorMin);
      foreach $l(@line_gaps){
        if(abs($median_value-$l)<=$radius){
          push(@line_gaps_new,$l)
        }
      }

      $num_lines_new=scalar(@line_gaps_new);
      if($num_lines_new==1){
        $exit_code=-1;
        last;
      }elsif($num_lines_new==2){
        $new_median_value=abs($line_gaps_new[0]+$line_gaps_new[1])/2;
        $new_median_value+=0.000001 if($new_median_value==0);
      }else{
        $new_median_value=$line_gaps_new[int(scalar(@line_gaps_new)/2)];
        $new_median_value+=0.000001 if($new_median_value==0);
      }
    }
    #print "DEBUG $exit_code $new_median_value $radius\n";
    if($exit_code==0){
      for($i=0;$i<$#lines_sorted;$i++){
        if($line_gaps[$i]>=$new_median_value-$radius && $line_gaps[$i]<=$new_median_value+$radius){
          print $lines_sorted[$i]," 1\n";
        }else{
          print $lines_sorted[$i]," 0\n";
        }
      }
    }else{
      foreach $l(@lines_sorted){
        print $l," 0\n";
      }
    }
  }
}


sub byGap
{
  my @fa=split(/\s+/,$a);
  my @fb=split(/\s+/,$b);
  return ($fa[1] <=> $fb[1]);
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
