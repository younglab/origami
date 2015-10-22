#!/bin/bash

## structure of script influenced by Diego's work

sed -e "s/\/[0-9];[0-9]//" $1 > left_reads.fq
sed -e "s/\/[0-9];[0-9]//" $2 > right_reads.fq

cutadapt -n 3 -m 17 --overlap 10 -a forward="ACGCGATATCTTATCTGACT" -a reverse="AGTCAGATAAGATATCGCGT" -o left_trim.fq --untrimmed-output left_untrim.fq -p right_trim.fq --untrimmed-paired-output right_untrim.fq left_reads.fq right_reads.fq
