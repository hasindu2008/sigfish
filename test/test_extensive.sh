#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

EVALUATE(){
    ex ./sigfish eval ${REF_PAF} ${MY_PAF} > $EVAL || die "Running the tool failed"
    MAPPED=$(grep -w "mapped_testset" $EVAL | head -1 | awk '{print $3}' | tr -d '(%)' )
    CORRECT=$(grep -w "correct" $EVAL | head -1 | awk '{print $3}' | tr -d '(%)')
    if (( $(echo "$MAPPED < $MAPPED_THRSH" | bc -l) )); then
        die "Test Failed: Mapped percentage ($MAPPED%) is less than $MAPPED_THRSH%"
    else
        echo "Test Passed: Mapped percentage $MAPPED%"
    fi

    if (( $(echo "$CORRECT < $CORRECT_THRSH" | bc -l) )); then
        die "Test Failed: Correct percentage ($CORRECT%) is less than $CORRECT_THRSH%"
    else
        echo "Test Passed: correct percentage $CORRECT%"
    fi

    echo "____________________________________________________"
    echo ""

}

THREADS=8

make

REF=test/nCoV-2019.reference.fasta
BLOW5=test/sp1_dna/batch0.blow5
REF_PAF=test/sp1_dna/batch0.minimap2.paf
MY_PAF=test/sp1_dna/test.paf
EVAL=test/sp1_dna/eval.txt
MAPPED_THRSH=100.0
CORRECT_THRSH=85.0

echo "DNA sp1"
ex  ./sigfish dtw ${REF} ${BLOW5} -t ${THREADS}  > ${MY_PAF} || die "Running the DTW failed"
EVALUATE

echo "DNA sp1 - from end"
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --from-end > ${MY_PAF} || die "Running DTW failed"
EVALUATE

REF=test/rnasequin_sequences_2.4.fa
BLOW5=test/sequin_rna/sequin_reads.blow5
REF_PAF=test/sequin_rna/rna.minimap2.paf
MY_PAF=test/sequin_rna/test.paf
EVAL=test/sequin_rna/eval.txt
MAPPED_THRSH=100.0
CORRECT_THRSH=75.0

echo "RNA sequin"
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna -q 500  -p -1 > ${MY_PAF}  || die "Running the DTW failed"
EVALUATE

echo "RNA sequin - full ref"
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --full-ref -q 500 -p -1 > ${MY_PAF}  || die "Running the DTW failed"
EVALUATE

echo "RNA sequin - from end"
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --from-end -q 500 > ${MY_PAF}  || die "Running the DTW failed"
EVALUATE

echo "RNA sequin - full ref from end"
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna --full-ref --from-end -q 500 > ${MY_PAF}  || die "Running the DTW failed"
EVALUATE


# echo "RNA invert"
# ./sigfish dtw ${REF} ${BLOW5} --rna --full-ref --invert -q 100000 -p0 -K10 --debug-break=yes > test/sequin_rna/res_rna.paf  || die "Running the DTW failed"
# ./sigfish dtw test/rnasequin_sequences_2.4.fa test/sequin_rna/sequin_reads.blow5 --rna --full-ref -q100000 -p0 -K10 --debug-break=yes> test/sequin_rna/res_rna2.paf  || die "Running the DTW failed"
# diff test/sequin_rna/res_rna.paf test/sequin_rna/res_rna2.paf || die "DTW results are different"

echo "*******************************************************"
echo "Tests passed"
echo "*******************************************************"

