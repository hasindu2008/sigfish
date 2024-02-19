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

handle_tests() {
	numfailed=$(wc -l < diff.txt)
	numcases=$(wc -l < ${ORIG})
	numres=$(wc -l < ${RES})
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt ${THRESH} ] && die "${1}: Validation failed"
	echo "Validation passed"
}

execute_test() {
	ORIG=$1
	RES=$2
	THRESH=$3
	diff -y --suppress-common-lines ${ORIG} ${RES} > diff.txt || handle_tests $testdir
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

# dtw

REF=test/nCoV-2019.reference.fasta
BLOW5=test/sp1_dna.blow5
REF_PAF=test/sp1_dna.minimap2.paf
MY_PAF=test/test.paf
EVAL=test/eval.txt
MAPPED_THRSH=100.0
CORRECT_THRSH=85.0

echo "DNA sp1"
ex  ./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} > ${MY_PAF} || die "Running the DTW failed"
EVALUATE

REF=test/rnasequin_sequences_2.4.fa
BLOW5=test/sequin_rna.blow5
REF_PAF=test/sequin_rna.minimap2.paf
MY_PAF=test/test.paf
EVAL=test/eval.txt
MAPPED_THRSH=100.0
CORRECT_THRSH=75.0

echo "RNA sequin"
ex ./sigfish dtw ${REF} ${BLOW5} -t ${THREADS} --rna -q 500  -p -1 > ${MY_PAF}  || die "Running the DTW failed"
EVALUATE

echo "RNA real jnn+dtw multi-thread"
ex ./sigfish real ${REF} ${BLOW5} -t 1 > ${MY_PAF} || die "Running the tool failed"
EVALUATE

echo "RNA real jnn+dtw multi-thread"
ex ./sigfish real ${REF} ${BLOW5} -t ${THREADS} > ${MY_PAF} || die "Running the tool failed"
EVALUATE

# realtime prefix

echo "RNA real prefix"
ex ./sigfish real test/sequin_rna.blow5 > test/prefix_real_rna.txt || die "Running the tool failed"
execute_test test/prefix_real_rna.txt test/data/prefix_real_rna.exp 5 || die "diff failed"

echo "DNA real prefix"
ex ./sigfish real test/sp1_dna.blow5 > test/prefix_real_dna.txt || die "Running the tool failed"
execute_test test/prefix_real_dna.txt test/data/prefix_real_dna.exp 5 || die "diff failed"

echo "*******************************************************"
echo "Tests passed"
echo "*******************************************************"

