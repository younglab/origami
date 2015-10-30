#!/bin/bash

OUTDIR=$1

macs2 callpeak -t $OUTDIR/mapped_reads.bam -n peaks --nomodel --shift 100 --outdir $OUTDIR
