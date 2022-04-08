#!/bin/bash

set -e


make
./sigfish dtw -g test/rnasequin_sequences_2.4.fa -s test/sequin_reads.blow5 --rna --full-ref \
    --invert -q100000 -p0 -K10 --debug-break=yes > test/res_rna.paf
./sigfish dtw -g test/rnasequin_sequences_2.4.fa -s test/sequin_reads.blow5 --rna --full-ref \
             -q100000 -p0 -K10 --debug-break=yes> test/res_rna2.paf
diff test/res_rna.paf test/res_rna2.paf
