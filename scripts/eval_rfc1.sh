#!/bin/bash

set -e

# tools locations
HARU_VENV=/mnt/d/Projects/HARU/RUscripts/env
MINIMAP2=minimap2

# test settings
RUN_MINIMAP2=${RUN_MINIMAP2:-true}
RUN_FPGA=${RUN_FPGA:-0}
TEST_CASE=${TEST_CASE:-1}
THREADS=8


# test data
TESTDATA_DIR=test/rfc1
RESULTS_DIR=test/rfc1/results
REF=${TESTDATA_DIR}/rfc1.fa

function form_mypaf_name () {
    if [ $RUN_FPGA -eq 1 ]; then
        MY_PAF="${MY_PAF}_fpga".paf
    else
        MY_PAF="${MY_PAF}".paf
    fi
}

function my_cmd () {
    $@
}

#${MINIMAP2} -cx map-ont ${REF} ${FASTQ} --secondary=no -t ${THREADS} > ${REF_PAF}
#./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --from-end > ${MY_PAF}
echo "mkdir -p ${TESTDATA_DIR}/run"
mkdir -p ${TESTDATA_DIR}/results

echo "Sourcing ${HARU_VENV}/bin/activate"
source ${HARU_VENV}/bin/activate

echo "TEST_CASE: ${TEST_CASE}"
if [ $TEST_CASE -eq 1 ]; then
    # batch0

    echo "=========================================================="
    echo "== Running test for batch0 no_rfc1"
    echo "=========================================================="
    BLOW5=${TESTDATA_DIR}/no_rfc1.blow5
    MY_PAF=${RESULTS_DIR}/test_no_rfc1
    form_mypaf_name
    REF_PAF=${TESTDATA_DIR}/no_rfc1_reference.paf
    my_cmd ./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} > ${MY_PAF}
    my_cmd uncalled pafstats -r ${REF_PAF} ${MY_PAF} > ${MY_PAF}.stats
    cat ${MY_PAF}.stats
    echo "Stats store in ${MY_PAF}.stats"

    echo "=========================================================="
    echo "== Running test for batch0 rfc1"
    echo "=========================================================="
    BLOW5=${TESTDATA_DIR}/rfc1.blow5
    MY_PAF=${RESULTS_DIR}/test_rfc1
    form_mypaf_name
    REF_PAF=${TESTDATA_DIR}/rfc1_reference.paf
    my_cmd ./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS}  > ${MY_PAF}
    my_cmd uncalled pafstats -r ${REF_PAF} ${MY_PAF} > ${MY_PAF}.stats
    cat ${MY_PAF}.stats
    
    echo "[INFO] Reference paf: ${REF_PAF}"
    echo "[INFO] My paf:        ${MY_PAF}"
    echo "[INFO] Stats:         ${MY_PAF}.stats"
elif [ $TEST_CASE -eq 2 ]; then
    # full sp1
    echo "Running test for full sp1 dataset"
    BLOW5=$2
    MY_PAF=${RESULTS_DIR}/test_sp1
    REF_PAF=${TESTDATA_DIR}/batch0.minimap2.paf
fi
