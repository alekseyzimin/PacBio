#!/bin/bash
CHR=$1
REF=$2
REF_NOISE=$3
SCF=$4
OVERHANG=1000000
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}


function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}
log "Working on chromosome $CHR in $CHR.dir"
mkdir -p $CHR.dir && \
ufasta extract -n chr$CHR $REF |awk '{print tolower($0)}' > $CHR.dir/ref$CHR.fa && \
ufasta extract -n chr$CHR $SCF >  $CHR.dir/qry$CHR.scaf.fa && \
(cd $CHR.dir && \
  mkdir -p  pass1 && \
  mkdir -p  pass2 && \
  if [ ! -s $CHR.all.polished.fa ];then
    splitScaffoldsAtNs.sh qry$CHR.scaf.fa 1  >  qry$CHR.fa.tmp && mv qry$CHR.fa.tmp qry$CHR.fa && \
    NUM_CONTIGS=`grep '^>' qry$CHR.fa | wc -l |awk '{print $1}'`
    if [ $NUM_CONTIGS = 1 ];then
      log "Single contig in chr$CHR, skipping pass1" && \
      mv qry$CHR.scaf.fa pass1/qry$CHR.scaf.fa.split.joined.fa
    else
      mkdir -p  pass1 && \
      mkdir -p  pass2 && \
      (cd pass1 && close_scaffold_gaps.sh -d asm -r ../qry$CHR.scaf.fa -q ../ref$CHR.fa -t 4 -m 4000 -o $OVERHANG 1>out.err 2>&1)
    fi
    if [ -s pass1/qry$CHR.scaf.fa.split.joined.fa ];then
        log "Running pass2 on chr$CHR" && \
        (cd pass2 && \
        ufasta one ../ref$CHR.fa |perl -ane 'BEGIN{$keep_on_ends=100000;}{if($F[0] =~ /^>/){print}else{$sub=length($F[0])-2*$keep_on_ends;if($sub>10000){substr($F[0],$keep_on_ends,$sub,"N"x$sub); print $F[0],"\n";}else{print}}}' | \
        splitScaffoldsAtNs.pl > ref$CHR.BEG_END.fa && \
        cat <(perl -e '{print ">chr'$CHR'\n"}') \
        <(grep -v '^>' ref$CHR.BEG_END.fa |head -n 1) \
        <(perl -e '{print "N"x1000}') \
        <(grep -v '^>' ../pass1/qry$CHR.scaf.fa.split.joined.fa ) \
        <(perl -e '{print "N"x1000}') \
        <(grep -v '^>' ref$CHR.BEG_END.fa |tail -n 1) > qry$CHR.scaf.padded.fa.tmp && mv qry$CHR.scaf.padded.fa.tmp qry$CHR.scaf.padded.fa && \
        splitScaffoldsAtNs.pl  < qry$CHR.scaf.padded.fa | ufasta one > asm$CHR.split.fa.tmp && mv asm$CHR.split.fa.tmp asm$CHR.split.fa && \
        grep '^>' --text asm$CHR.split.fa | perl -ane '{($rn,$coord)=split(/\./,substr($F[0],1));$h{$rn}.=substr($F[0],1)." ";}END{foreach $r(keys %h){@f=split(/\s+/,$h{$r}); for ($i=0;$i<$#f;$i++){print $f[$i]," ",$f[$i+1],"\n"}}}' > asm$CHR.split.fa.valid_join_pairs.tmp && mv asm$CHR.split.fa.valid_join_pairs.tmp asm$CHR.split.fa.valid_join_pairs && \
        nucmer --batch $OVERHANG -p aln -l 21 -c 200 -t 4 <(ufasta one asm$CHR.split.fa |perl -ane 'BEGIN{$keep_on_ends='$OVERHANG';}{if($F[0] =~ /^>/){print}else{$sub=length($F[0])-2*$keep_on_ends;if($sub>10000){substr($F[0],$keep_on_ends,$sub,"N"x$sub); print $F[0],"\n";}else{print}}}') ../ref$CHR.fa && \
        delta-filter -1 aln.delta > aln.1delta.tmp && mv aln.1delta.tmp aln.1delta && \
        show-coords -lcHq aln.1delta | extract_merges.pl ../ref$CHR.fa 2500 $OVERHANG asm asm$CHR.split.fa.valid_join_pairs > $CHR.links.tmp && mv $CHR.links.tmp $CHR.links && \
        merge_contigs.pl asm$CHR.split.fa < $CHR.links | create_merged_sequences.pl asm$CHR.split.fa  $CHR.links > $CHR.scaffolds.fa.tmp && mv $CHR.scaffolds.fa.tmp $CHR.scaffolds.fa && \
        recover_scaffolds.pl < $CHR.scaffolds.fa |ufasta format > $CHR.joined.fa.tmp && mv $CHR.joined.fa.tmp ../$CHR.all.polished.fa )
      else
        log "Warning! Unable to run pass2 on chr$CHR" && \
        cp qry$CHR.scaf.fa $CHR.all.polished.fa.tmp && mv $CHR.all.polished.fa.tmp $CHR.all.polished.fa
      fi
  fi && \
  log "Success!!! Output sequences in $CHR.all.polished.fa"
)
