#!/bin/bash

source ~/dsday/origami/bin/dispatch.sh

OUTDIR=$1
PARALLEL=$2

## structure of script influenced by Diego's work

dispatch "sed -e 's/\/[0-9];[0-9]//' $3 > $OUTDIR/left_reads.fq"
dispatch "sed -e 's/\/[0-9];[0-9]//' $4 > $OUTDIR/right_reads.fq"

wait

dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/l_t1.fq --untrimmed-output $OUTDIR/l_nt1.fq -p $OUTDIR/r_t1.fq --untrimmed-paired-output $OUTDIR/r_nt1.fq $OUTDIR/left_reads.fq $OUTDIR/right_reads.fq
wait
dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/r_t2.fq --untrimmed-output $OUTDIR/r_nt2.fq -p $OUTDIR/l_t2.fq --untrimmed-paired-output $OUTDIR/l_nt2.fq $OUTDIR/r_nt1.fq $OUTDIR/l_nt1.fq
dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/r_t3.fq --untrimmed-output $OUTDIR/r_nt3.fq -p $OUTDIR/l_t3.fq --untrimmed-paired-output $OUTDIR/l_nt3.fq $OUTDIR/r_t1.fq $OUTDIR/l_t1.fq

wait

cat $OUTDIR/l_t3.fq $OUTDIR/l_nt3.fq $OUTDIR/l_t2.fq > $OUTDIR/left_kept.fq
cat $OUTDIR/r_t3.fq $OUTDIR/r_nt3.fq $OUTDIR/r_t2.fq > $OUTDIR/right_kept.fq

cat $OUTDIR/l_nt2.fq > $OUTDIR/left_untrimmed.fq
cat $OUTDIR/r_nt2.fq > $OUTDIR/right_untrimmed.fq

### Cleanup
rm $OUTDIR/l_*.fq $OUTDIR/r_*.fq
rm $OUTDIR/left_reads.fq $OUTDIR/right_reads.fq
