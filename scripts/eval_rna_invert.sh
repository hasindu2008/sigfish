#!/bin/bash

set -e


make
./sigfish dtw test/rnasequin_sequences_2.4.fa test/sequin_rna/sequin_reads.blow5 --rna --full-ref \
    --invert -q100000 -p0 -K10 --debug-break=yes > test/sequin_rna/res_rna.paf
./sigfish dtw test/rnasequin_sequences_2.4.fa test/sequin_rna/sequin_reads.blow5 --rna --full-ref \
             -q100000 -p0 -K10 --debug-break=yes> test/sequin_rna/res_rna2.paf
diff test/sequin_rna/res_rna.paf test/sequin_rna/res_rna2.paf
