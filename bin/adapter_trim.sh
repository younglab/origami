#!/bin/bash

source ~/dsday/origami/bin/dispatch.sh

OUTDIR=$1
PARALLEL=$2
SPLITNUM=$3

## structure of script influenced by Diego's work

if [ $PARALLEL = "on" ]
then
  # need to generalize this
  dispatch "bzcat $4 | sed -e 's/\/[0-9];[0-9]//' > $OUTDIR/tmp/left_reads.fq"
  dispatch "bzcat $5 | sed -e 's/\/[0-9];[0-9]//' > $OUTDIR/tmp/right_reads.fq"
  
  wait 
  dispatch "split -l $SPLITNUM $OUTDIR/tmp/left_reads.fq $OUTDIR/tmp/leftreads"
  dispatch "split -l $SPLITNUM $OUTDIR/tmp/right_reads.fq $OUTDIR/tmp/rightreads"
  
  wait
  
  #rm $OUTDIR/tmp/left_reads.fq $OUTDIR/tmp/right_reads.fq
  dispatch "bzip2 $OUTDIR/tmp/left_reads.fq"
  dispatch "bzip2 $OUTDIR/tmp/right_reads.fq"

  wait
  
  ## One assumption here is that split names the files in the same linear order -- maybe this should be done differently?
  LEFTREADS=($(ls $OUTDIR/tmp/leftreads*))
  RIGHTREADS=($(ls $OUTDIR/tmp/rightreads*))
  for ((i=0;i<${#LEFTREADS[@]};++i)); do
    dispatch cutadapt -f fastq -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/tmp/l_t1_$i.fq --untrimmed-output $OUTDIR/tmp/l_nt1_$i.fq -p $OUTDIR/tmp/r_t1_$i.fq --untrimmed-paired-output $OUTDIR/tmp/r_nt1_$i.fq ${LEFTREADS[$i]} ${RIGHTREADS[$i]}
  done
  wait
  
  for ((i=0;i<${#LEFTREADS[@]};++i)); do
    dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/tmp/r_t2_$i.fq --untrimmed-output $OUTDIR/tmp/r_nt2_$i.fq -p $OUTDIR/tmp/l_t2_$i.fq --untrimmed-paired-output $OUTDIR/tmp/l_nt2_$i.fq $OUTDIR/tmp/r_nt1_$i.fq $OUTDIR/tmp/l_nt1_$i.fq
    dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/tmp/r_t3_$i.fq --untrimmed-output $OUTDIR/tmp/r_nt3_$i.fq -p $OUTDIR/tmp/l_t3_$i.fq --untrimmed-paired-output $OUTDIR/tmp/l_nt3_$i.fq $OUTDIR/tmp/r_t1_$i.fq $OUTDIR/tmp/l_t1_$i.fq
  done

  wait
  
  dispatch "cat $OUTDIR/tmp/l_t3*.fq $OUTDIR/tmp/l_nt3*.fq $OUTDIR/tmp/l_t2*.fq > $OUTDIR/tmp/left_kept.fq"
  dispatch "cat $OUTDIR/tmp/r_t3*.fq $OUTDIR/tmp/r_nt3*.fq $OUTDIR/tmp/r_t2*.fq > $OUTDIR/tmp/right_kept.fq"

  dispatch "cat $OUTDIR/tmp/l_nt2*.fq > $OUTDIR/tmp/left_untrimmed.fq"
  dispatch "cat $OUTDIR/tmp/r_nt2*.fq > $OUTDIR/tmp/right_untrimmed.fq"

  wait
  
  rm $OUTDIR/tmp/leftreads* $OUTDIR/tmp/rightreads*

else

  dispatch "sed -e 's/\/[0-9];[0-9]//' $4 > $OUTDIR/tmp/left_reads.fq"
  dispatch "sed -e 's/\/[0-9];[0-9]//' $5 > $OUTDIR/tmp/right_reads.fq"

  wait

  dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/tmp/l_t1.fq --untrimmed-output $OUTDIR/tmp/l_nt1.fq -p $OUTDIR/tmp/r_t1.fq --untrimmed-paired-output $OUTDIR/tmp/r_nt1.fq $OUTDIR/tmp/left_reads.fq $OUTDIR/tmp/right_reads.fq
  wait
  dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/tmp/r_t2.fq --untrimmed-output $OUTDIR/tmp/r_nt2.fq -p $OUTDIR/tmp/l_t2.fq --untrimmed-paired-output $OUTDIR/tmp/l_nt2.fq $OUTDIR/tmp/r_nt1.fq $OUTDIR/tmp/l_nt1.fq
  dispatch cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o $OUTDIR/tmp/r_t3.fq --untrimmed-output $OUTDIR/tmp/r_nt3.fq -p $OUTDIR/tmp/l_t3.fq --untrimmed-paired-output $OUTDIR/tmp/l_nt3.fq $OUTDIR/tmp/r_t1.fq $OUTDIR/tmp/l_t1.fq

  wait

  dispatch "cat $OUTDIR/tmp/l_t3.fq $OUTDIR/tmp/l_nt3.fq $OUTDIR/tmp/l_t2.fq > $OUTDIR/tmp/left_kept.fq"
  dispatch "cat $OUTDIR/tmp/r_t3.fq $OUTDIR/tmp/r_nt3.fq $OUTDIR/tmp/r_t2.fq > $OUTDIR/tmp/right_kept.fq"

  dispatch "cat $OUTDIR/tmp/l_nt2.fq > $OUTDIR/tmp/left_untrimmed.fq"
  dispatch "cat $OUTDIR/tmp/r_nt2.fq > $OUTDIR/tmp/right_untrimmed.fq"

  wait
  
  rm $OUTDIR/tmp/left_reads.fq $OUTDIR/tmp/right_reads.fq
  
fi

### Cleanup
rm $OUTDIR/tmp/l_*.fq $OUTDIR/tmp/r_*.fq

dispatch bzip2 $OUTDIR/tmp/left_untrimmed.fq
dispatch bzip2 $OUTDIR/tmp/right_untrimmed.fq

wait

#rm -f $OUTDIR/tmp/left_untrimmed.fq $OUTDIR/tmp/right_untrimmed.fq

