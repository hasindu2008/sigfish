#!/bin/bash

set -e

HARU_VENV=~/haru/
REF=test/rnasequin_sequences_2.4.fa
BLOW5=test/sequin_rna/sequin_reads.blow5
THREADS=8
REF_PAF=test/sequin_rna/rna.minimap2.paf
MY_PAF=test/sequin_rna/test.paf

make
#./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --from-end -q 500 > ${MY_PAF}
#./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --full-ref --from-end -q 500 > ${MY_PAF}
#./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --full-ref -q 500  > ${MY_PAF}
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna -q 500  -p -1 > ${MY_PAF}
#./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --full-ref -q 500 -p -1 > ${MY_PAF}

source ${HARU_VENV}/bin/activate
uncalled pafstats -r ${REF_PAF} ${MY_PAF} -a > test/sequin_rna/err.paf
./sigfish eval ${REF_PAF} ${MY_PAF}