#!/bin/bash

set -e


if [ -z "$1" ] || [ -z "$2"]; then
    echo "Usage: $0 ref_id read_id"
    exit 1
fi

samtools faidx test/rnasequin_sequences_2.4.fa $1 > a.fa
./sigfish dtw -g a.fa -s test/sequin_reads.blow5 -t 16 --rna --full-ref -q 500 -p -1 > a.paf
grep $2 a.paf



