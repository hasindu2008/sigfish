#!/bin/bash

set -e

HARU_VENV=~/haru/
REF=test/rnasequin_sequences_2.4.fa
BLOW5=test/sequin_reads.blow5
THREADS=8
REF_PAF=test/rna.minimap2.paf
MY_PAF=test/test.paf

make
#./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --rna --from-end -q 500 > ${MY_PAF}
#./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --rna --full-ref --from-end -q 500 > ${MY_PAF}
#./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --rna --full-ref -q 500  > ${MY_PAF}
./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --rna -q 500  -p -1 > ${MY_PAF}
#./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --rna --full-ref -q 500 -p -1 > ${MY_PAF}

./sigfish eval ${REF_PAF} ${MY_PAF}
source ${HARU_VENV}/bin/activate
uncalled pafstats -r ${REF_PAF} ${MY_PAF}
uncalled pafstats -r ${REF_PAF} ${MY_PAF} -a > err.paf
./sigfish eval ${REF_PAF} ${MY_PAF}