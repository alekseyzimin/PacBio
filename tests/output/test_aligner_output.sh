#! /bin/sh

SR=test_super_reads.fa
PB=test_pacbio.fa
UL=test_unitigs_lengths

RETURN=0
compare() {
    if diff -q $1 $1_expected; then
        echo "[OK  ] $1"
    else
        wdiff $1 $1_expected
        echo "[FAIL] $1"
        RETURN=1
    fi
}

../../src_jf_aligner/jf_aligner -s 10k -m 17 -r $SR -p $PB --details details_normal --coords coords_normal || RETURN=1
compare coords_normal
compare details_normal

../../src_jf_aligner/jf_aligner -s 10k -m 17 -r $SR -p $PB -f -l $UL -k 65 --details details_forward --coords coords_forward || RETURN=1
compare coords_forward
compare details_forward

exit $RETURN
