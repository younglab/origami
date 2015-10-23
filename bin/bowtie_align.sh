#!/bin/bash

source ~/dsday/origami/bin/dispatch.sh

OUTDIR=$1
PARALLEL=$2

dispatch "bowtie -n 1 -m 1 -p 8 --sam /nfs/genomes/human_gp_feb_09/bowtie/hg19 $OUTDIR/tmp/left_kept.fq > $OUTDIR/tmp/left_kept.sam"
dispatch "bowtie -n 1 -m 1 -p 8 --sam /nfs/genomes/human_gp_feb_09/bowtie/hg19 $OUTDIR/tmp/right_kept.fq > $OUTDIR/tmp/right_kept.sam"

wait

dispatch "samtools view -Sb $OUTDIR/tmp/left_kept.sam > $OUTDIR/tmp/left_kept.bam"
dispatch "samtools view -Sb $OUTDIR/tmp/right_kept.sam > $OUTDIR/tmp/right_kept.bam"

wait

rm $OUTDIR/tmp/left_kept.sam $OUTDIR/tmp/right_kept.sam
mv $OUTDIR/tmp/left_kept.bam $OUTDIR/left_kept.bam
mv $OUTDIR/tmp/right_kept.bam $OUTDIR/right_kept.bam
