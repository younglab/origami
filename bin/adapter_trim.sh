#!/bin/bash

OUTDIR=$1

## structure of script influenced by Diego's work

sed -e "s/\/[0-9];[0-9]//" $2 > $OUTDIR/left_reads.fq
sed -e "s/\/[0-9];[0-9]//" $3 > $OUTDIR/right_reads.fq

cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/l_t1.fq --untrimmed-output $OUTDIR/l_nt1.fq -p $OUTDIR/r_t1.fq --untrimmed-paired-output $OURDIR/r_nt1.fq $OURDIR/left_reads.fq $OUTDIR/right_reads.fq
cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/r_t2.fq --untrimmed-output $OUTDIR/r_nt2.fq -p $OUTDIR/l_t2.fq --untrimmed-paired-output $OURDIR/l_nt2.fq $OUTDIR/r_nt1.fq $OUTDIR/l_nt1.fq
cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/r_t3.fq --untrimmed-output $OUTDIR/r_nt3.fq -p $OUTDIR/l_t3.fq --untrimmed-paired-output $OURDIR/l_nt3.fq $OUTDIR/r_t1.fq $OUTDIR/l_t1.fq

cat $OUTDIR/l_t3.fq $OUTDIR/l_nt3.fq $OURDIR/l_t2.fq > $OUTDIR/left_kept.fq
cat $OUTDIR/r_t3.fq $OUTDIR/r_nt3.fq $OURDIR/r_t2.fq > $OUTDIR/right_kept.fq

cat $OUTDIR/l_nt2.fq > $OUTDIR/left_untrimmed.fq
cat $OUTDIR/r_nt2.fq > $OUTDIR/right_untrimmed.fq

