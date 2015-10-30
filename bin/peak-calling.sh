#!/bin/bash

OUTDIR=$1

bedtools bamtobed -i $OUTDIR/mapped_reads.bam > $OUTDIR/tmp/reads.bed
macs2 callpeak -t $OUTDIR/tmp/reads.bed -n peaks --nomodel --shift 100 --outdir $OUTDIR
rm $OUTDIR/tmp/reads.bed
