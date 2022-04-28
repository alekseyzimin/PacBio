#!/bin/bash
chrM_REF=$1
ASM=$2
nucmer -c 200 -t 32 -p chrM $chrM_REF $ASM
ASM_ctg=`delta-filter -r chrM.delta |show-coords -lcHr /dev/stdin | perl -ane '{$ctg=$F[-1];if($F[3]<$F[4]){$c1=$F[3];$c2=$F[4];}else{$c1=$F[4];$c2=$F[3];}END{print $ctg,"\n"}}'`
ASM_beg=`delta-filter -r chrM.delta |show-coords -lcHr /dev/stdin | perl -ane '{$ctg=$F[-1];if($F[3]<$F[4]){$c1=$F[3];$c2=$F[4];}else{$c1=$F[4];$c2=$F[3];}END{print $c1,"\n"}}'`
ASM_end=`delta-filter -r chrM.delta |show-coords -lcHr /dev/stdin | perl -ane '{$ctg=$F[-1];if($F[3]<$F[4]){$c1=$F[3];$c2=$F[4];}else{$c1=$F[4];$c2=$F[3];}END{print $c2,"\n"}}'`
ufasta extract -n $ASM_ctg $ASM |ufasta one |grep -v '^>' |awk '{print ">chrM\n"substr($1,'$ASM_beg',('$ASM_end'-'$ASM_beg'+1))}' > chrM.new.fa
