#!/bin/bash

source ~/dsday/origami/bin/dispatch.sh

OUTDIR=$1
PARALLEL=$2

# Split into 1 million reads per file (1 FASTQ read = 4 lines)
dispatch "split -l 4000000 $OUTDIR/tmp/left_kept.fq $OUTDIR/tmp/leftkept"
dispatch "split -l 4000000 $OUTDIR/tmp/right_kept.fq $OUTDIR/tmp/rightkept"

wait

for FILE in $OUTDIR/tmp/leftkept*
do
	dispatch "bowtie -n 1 -m 1 -p 6 --sam /nfs/genomes/human_gp_feb_09/bowtie/hg19 $FILE > $FILE.sam; samtools view -Sb $FILE.sam > $FILE.bam; rm $FILE.sam"
done

for FILE in $OUTDIR/tmp/rightkept*
do
	dispatch "bowtie -n 1 -m 1 -p 6 --sam /nfs/genomes/human_gp_feb_09/bowtie/hg19 $FILE > $FILE.sam; samtools view -Sb $FILE.sam > $FILE.bam; rm $FILE.sam"
done

wait

dispatch "samtools merge $OUTDIR/tmp/left_kept.bam $OUTDIR/tmp/leftkept*.bam"
dispatch "samtools merge $OUTDIR/tmp/right_kept.bam $OUTDIR/tmp/rightkept*.bam"

wait

dispatch "rm $OUTDIR/tmp/leftkept* $OUTDIR/tmp/rightkept*"
wait

dispatch "samtools sort -Obam -Tlefttmp -n $OUTDIR/tmp/left_kept.bam > $OUTDIR/tmp/left_kept.sorted.bam"
dispatch "samtools sort -Obam -Trighttmp -n $OUTDIR/tmp/right_kept.bam > $OUTDIR/tmp/right_kept.sorted.bam"

wait

dispatch "~/dsday/origami/bin/mapped-reads-merge $OUTDIR/tmp/left_kept.sorted.bam $OUTDIR/tmp/right_kept.sorted.bam $OUTDIR/tmp/mapped_reads.bam"

wait

dispatch "samtools sort -Obam -Ttmp $OUTDIR/tmp/mapped_reads.bam > $OUTDIR/mapped_reads.bam"

wait

rm $OUTDIR/tmp/left_kept.sam $OUTDIR/tmp/right_kept.sam
rm $OUTDIR/tmp/left_kept.sorted.bam $OUTDIR/tmp/right_kept.sorted.bam
rm $OUTDIR/tmp/mapped_reads.bam

#mv $OUTDIR/tmp/left_kept.bam $OUTDIR/left_kept.bam
#mv $OUTDIR/tmp/right_kept.bam $OUTDIR/right_kept.bam
