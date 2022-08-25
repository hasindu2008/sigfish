#!/bin/bash

set -e


make
./sigfish dtw test/rnasequin_sequences_2.4.fa test/rna_sequin/sequin_reads.blow5 --rna --full-ref \
    --invert -q100000 -p0 -K10 --debug-break=yes > test/rna_sequin/res_rna.paf
./sigfish dtw test/rnasequin_sequences_2.4.fa test/rna_sequin/sequin_reads.blow5 --rna --full-ref \
             -q100000 -p0 -K10 --debug-break=yes> test/rna_sequin/res_rna2.paf
diff test/rna_sequin/res_rna.paf test/rna_sequin/res_rna2.paf
