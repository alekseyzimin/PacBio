#!/usr/bin/env perl
$rn="";
$max_len=65500;
$shooting_index=0;
if($ARGV[0] eq ""){
$suffix="super-read";
}else{
$suffix=$ARGV[0];
}
print "$rn\n";
while($line=<STDIN>){
    if($line =~ /^>/){
	if(not($rn eq "")){
	    $l=length($seq);
	    if($l<$max_len){
		print "$rn\n$seq\n";
	    }else{
		$num_pieces=int($l/$max_len)+1;
		$max_len_local=int($l/$num_pieces)+1;
		$offset=int(($max_len_local)/2)+1;
	    	for($i=0;$i<$num_pieces*2-1;$i++){
			print "$rn.",$i*$offset,".",$i*$offset+$max_len_local,"\n",substr($seq,$i*$offset,$max_len_local),"\n";
		}
	    }
	}
	chomp($line); 
	@l=split(/\s+/,$line);
        $rn=$l[0].":".$suffix;
	$seq="";
    }else{
	chomp($line); 
	$seq.=$line;
    }
}
$l=length($seq);
            if($l<$max_len){
                print "$rn\n$seq\n";
            }else{
                $max_len_local=int(($l-10000)/$l*$max_len);
                $offset=int(($max_len_local-1)/2);
                    for($i=0;$i<$l;$i+=$offset){
                        print "$rn.$max_len_local.$i\n",substr($seq,$i,$max_len_local),"\n";
                    }   
            }

