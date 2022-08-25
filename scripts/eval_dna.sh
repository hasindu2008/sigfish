#!/bin/bash

set -e

HARU_VENV=~/haru/
MINIMAP2=minimap2
FASTQ=test/sp1_dna/batch0.fastq
REF=test/nCoV-2019.reference.fasta
BLOW5=test/sp1_dna/batch0.blow5
THREADS=8
REF_PAF=test/sp1_dna/batch0.minimap2.paf
MY_PAF=test/sp1_dna/test.paf

make
#${MINIMAP2} -cx map-ont ${REF} ${FASTQ} --secondary=no -t ${THREADS} > ${REF_PAF}
#./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --from-end > ${MY_PAF}
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS}  > ${MY_PAF}

source ${HARU_VENV}/bin/activate
uncalled pafstats -r ${REF_PAF} ${MY_PAF}
./sigfish eval ${REF_PAF} ${MY_PAF}