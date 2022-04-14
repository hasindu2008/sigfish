#!/bin/bash

set -e

scripts/eval_rna.sh
grep "rf:Z:tp" err.paf  | awk '{print $1}' > correct.list
grep "rf:Z:fp" err.paf  | awk '{print $1}' > wrong.list

samtools faidx test/sequin_reads.fastq
awk '{print $1":"$2-50"-"$2}' test/sequin_reads.fastq.fai > reads.list
#awk '{print $1":1-50"}' test/sequin_reads.fastq.fai > reads.list
samtools faidx test/sequin_reads.fastq -r reads.list > reads_50.fastq
RNAfold reads_50.fastq | awk '{if(NR%3==1){printf $1"\t"} else if (NR%3==0){print $NF} }' |  tr -d '>()' > reads_50.fold.mnn

rm -f correct.fold.mnn
cat correct.list | while read p
do
    grep $p reads_50.fold.mnn >> correct.fold.mnn
done

rm -f wrong.fold.mnn
cat wrong.list | while read p
do
    grep $p reads_50.fold.mnn >> wrong.fold.mnn
done



