#!/bin/bash

OUTDIR=$1

bowtie -n 1 -m 1 -p 8 --sam /nfs/genomes/human_gp_feb_09/bowtie/hg19 $OUTDIR/tmp/left_kept.fq > $OUTDIR/tmp/left_kept.sam
bowtie -n 1 -m 1 -p 8 --sam /nfs/genomes/human_gp_feb_09/bowtie/hg19 $OUTDIR/tmp/right_kept.fq > $OUTDIR/tmp/right_kept.sam

samtools view -Sb $OUTDIR/tmp/left_kept.sam > $OUTDIR/tmp/left_kept.bam
samtools view -Sb $OUTDIR/tmp/right_kept.sam > $OUTDIR/tmp/right_kept.bam

rm $OUTDIR/tmp/left_kept.sam $OUTDIR/tmp/right_kept.sam
