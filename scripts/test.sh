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
BLOW5=test/sp1_dna.blow5
REF_PAF=test/sp1_dna.minimap2.paf
MY_PAF=test/test.paf
EVAL=test/eval.txt
MAPPED_THRSH=100.0
CORRECT_THRSH=85.0

echo "DNA sp1"
ex  ./sigfish dtw ${REF} ${BLOW5} -t ${THREADS}  > ${MY_PAF} || die "Running the DTW failed"
EVALUATE

REF=test/rnasequin_sequences_2.4.fa
BLOW5=test/sequin_rna.blow5
REF_PAF=test/sequin_rna.minimap2.paf
MY_PAF=test/test.paf
EVAL=test/eval.txt
MAPPED_THRSH=100.0
CORRECT_THRSH=75.0

echo "RNA sequin"
./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna -q 500  -p -1 > ${MY_PAF}  || die "Running the DTW failed"
EVALUATE


echo "*******************************************************"
echo "Tests passed"
echo "*******************************************************"

