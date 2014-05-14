#!/usr/bin/env perl
#build mega-reads for a single pb read
@mega_reads=();
@mega_reads_last_sr=();
@mega_reads_last_pos=();
@starts=();
@ends=();
@sr=();
@pb_cover;
while($line=<STDIN>){
    chomp($line);
    @f=split(/\s+/,$line);
    if(scalar(@pb_cover)==0){#record coverage of pacbio reads to distinguish between spurs and gaps in coverage
	for($i=0;$i<=$f[10];$i++){
		push(@pb_cover,0);
	}
    }
   
    push(@istarts,$f[0]-$f[2]);
    push(@iends,$f[1]+$f[10]-$f[3]); 
    push(@starts,$f[0]);
    push(@ends,$f[1]);
    push(@scores,$f[8]+$f[1]);#we include the last position and number of kmers MUST SUBTRACT THE PREVIOUS LAST COORDINATE
    push(@sr,$f[12]);
    push(@lines,$line);
    push(@contained,0);
    push(@used,0);
    push(@spur,0);
    for($i=$f[0];$i<=$f[1];$i++){
	$pb_cover[$i]++;
    }
}

#first we pre-process the data by removing "spurs"; everything that does not overlap on the right is a spur, except for the last super-read
$spur_found=1;
while($spur_found){
    $spur_found=0;
    for($i=0;$i<$#sr;$i++){
	next if($spur[$i]);
        next if($pb_cover[$ends[$i]]<=1); #not a spur if nothing covering the end according to PB matches
	for($j=$i+1;$j<=$#sr;$j++){
	    next if($spur[$j]);
	    next if($starts[$j]>$ends[$i]);
	    $t_ext=overlap_ext($sr[$i],$sr[$j]);
	    last if($t_ext>0);
	}
	if($j>$#sr){ #extension does not continue
	    print "Found spur $i $pb_cover[$ends[$i]] $starts[$i] $ends[$i] $sr[$i]\n";
	    $spur[$i]=1;
	    $spur_found=1;

	    for(my $k=$starts[$i];$k<=$ends[$i];$k++){
        	$pb_cover[$k]--;
    		}
	}
    }
}


#find the first non-spur super-read 
for($i=0;$i<$#sr;$i++){
    last if(not($spur[$i]));
}
exit if($i>$#sr);
#create first mega-read
push(@mega_reads,$sr[$i]);
push(@mega_reads_first_pos,$istarts[$i]);
push(@mega_reads_last_sr,$sr[$i]);
push(@mega_reads_last_pos,$ends[$i]);
push(@mega_reads_last_ipos,$iends[$i]);
push(@mega_reads_last_index,$i);
push(@mega_reads_num_sr,1);
$contained[$i]=1;
$used[$i]=1;
$ext_flag=1;
$mri=0;

#main mega-read extension loop -- loop until all mega-reads are created
while(1){
    print "Extending mega-read $mega_reads[$mri]\n";
    $ext_flag=1;
    while($ext_flag){#if there was extension
	$ext_flag=0;
#here we look for the best extension for current mega-read, must overlap and have the highest score=amount of extension+k-mer coverage, only consider overlapping extensions
	@extensions=();
	$max_ext=-1000000000;
	$max_ext_index=-1;
	$max_overlap=0;
	print "Looking for extensions for $mega_reads_last_pos[$mri] $mega_reads_last_sr[$mri]\n";
	for($i=$mega_reads_last_index[$mri];$i<=$#sr;$i++){
	    next if($starts[$i]>$mega_reads_last_pos[$mri]);#if pb indicates no overlap
            next if($contained[$i]);#if already declared contained
	    next if($spur[$i]);#if spur
	    $overlap=overlap_ext($mega_reads_last_sr[$mri],$sr[$i]);
	    if($overlap>0){
		$ext_score=$scores[$i]-$mega_reads_last_pos[$mri];#here we score the continuation
		print "New extension $i $ext_score $starts[$i] $ends[$i]  $sr[$i] $mega_reads_last_pos[$mri] $scores[$i]\n"; 
		push(@extensions,$i);
		if($ext_score>=$max_ext){
		    $max_ext=$ext_score;
		    $max_ext_index=$i;
		    $max_overlap=$overlap;
		}
	    }
	}
	print "Candidate max extension $max_ext\nindex $max_ext_index\nsr $sr[$max_ext_index]\n";
	if($max_ext_index>-1){
	    $ext_flag=1;
	    $used[$max_ext_index]=1;
#now eliminate all containees
	    @t=split(/_/,$sr[$max_ext_index]);
	    print "Old mega read $mega_reads[$mri]\n";
	    print "Extension     $sr[$max_ext_index]\n";
	    $mega_reads[$mri].="_".join("_",@t[$#t-$max_overlap+1..$#t]);
	    $mega_reads_num_sr[$mri]++;
	    print "Extended mega read $mega_reads[$mri]\n";
	    print "Checking containees\n";
	    $mega_reads_last_ipos[$mri]=$iends[$max_ext_index];
	    $mega_reads_last_pos[$mri]=$ends[$max_ext_index];
	    $mega_reads_last_sr[$mri]=$sr[$max_ext_index];
	    $mega_reads_last_index[$mri]=$max_ext_index;
	    for($i=0;$i<$#extensions;$i++){
		print "Checking containee $i $sr[$extensions[$i]]\n";
		$contained[$extensions[$i]]=1 if(index($mega_reads[$mri],$sr[$extensions[$i]])>-1);
		if(index($mega_reads[$mri],$sr[$extensions[$i]])>-1){
		    $contained[$extensions[$i]]=1;
		    print "Found containee $extensions[$i] $sr[$extensions[$i]]\n";
		}
	    }
	}
    }
#now we see if there is anything left to try
    for($i=0;$i<=$#sr;$i++){
	if(not($contained[$i]) && not($spur[$i]) && not($used[$i])){
	    print "New mega-read starting $sr[$i]\n";
	    $mri++;
	    $contained[$i]=1;
	    $used[$i]=1;
	    push(@mega_reads_first_pos,$istarts[$i]);
	    push(@mega_reads,$sr[$i]);
	    push(@mega_reads_last_ipos,$iends[$i]);
	    push(@mega_reads_last_sr,$sr[$i]);
	    push(@mega_reads_last_pos,$ends[$i]);
	    push(@mega_reads_last_index,$i);
	    push(@mega_reads_num_sr,1);
	    last;
	}
    }
    last if($i==$#sr+1);
}

#now we print all non-contained lines
for($i=0;$i<=$#lines;$i++){
    if($used[$i]){
	print "used $lines[$i]\n";
    }elsif($contained[$i]){
	print "contained $lines[$i]\n";
    }else{
	print "unknown $lines[$i]\n";
    }
}

for($i=0;$i<=$#mega_reads;$i++){
    print STDERR "$mega_reads_num_sr[$i] $mega_reads_first_pos[$i] $mega_reads_last_ipos[$i] ",$mega_reads_last_ipos[$i]-$mega_reads_first_pos[$i]," $mega_reads[$i]\n";
}

sub overlap_ext{
    my @f1=split(/_/,$_[1]);
    my $extension=0;
    for(my $j=$#f1;$j>=0;$j--){
        $tmpstr=join("_",@f1[0..$j]);
        if(index($_[0],$tmpstr)>0 && index($_[0],$tmpstr)==length($_[0])-length($tmpstr)){
            $extension=$#f1-$j;;
            last;
        }
    }
    return($extension);
}

