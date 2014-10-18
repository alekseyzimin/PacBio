#!/bin/bash
KMER=$2;
COORDS=$1;
awk 'BEGIN{flag=0}{
	if($0 ~ /^>/){
		flag=0;
		pb=substr($1,2);
	}else{
		flag++;
	}; 
	if(flag>1 && last_mr!=$8){
		l=split(last_mr,a,"_");
		split($8,b,"_");
		k1=int(substr(a[l],1,length(a[l])-1));
		k2=int(substr(b[1],1,length(b[1])-1));
		if(k1<k2){
			print pb" "$1-$3-last_coord" "k1" "k2
		}else 
			if(k1>k2){
				print pb" "$1-$3-last_coord" "k2" "k1
			}
	};
	last_mr=$8;
	last_coord=$2+$5-$4;
}' ${COORDS}.all.txt |sort -nk3 -k4n -S 10%|uniq -D -f 2 
