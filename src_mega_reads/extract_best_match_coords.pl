#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl
  
while($line=<STDIN>){
        chomp($line);
        $line=~s/^\s+//;
        @f=split(/\s+/,$line);
        $score=int($f[7]*$f[9]/100);
        if($score>$scores{$f[-1]}){
                $scores{$f[-1]}=$score;
                $lines{$f[-1]}=$line;
        }
}

foreach $v(keys %lines){
print $lines{$v},"\n";
}
