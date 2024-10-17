#! /bin/bash
# Test the bias when generating with an MC template and fitting with the parametric model

MCCARD=$1
DATACARD=$2
TAG=$3
TOYS=$4

if [[ "${MCCARD}" == "" ]] || [[ "${DATACARD}" == "" ]]; then
    echo "No input cards provided!"
    exit
fi
if [[ "${TOYS}" == "" ]]; then
    TOYS="500"
fi

${CMSSW_BASE}/src/ZLFV_fits/tests/bemu_gen_fit_test.sh --card_1 ${MCCARD} --card_2 ${DATACARD} --name ${TAG} -g 500 -t ${TOYS} -r 50 --skipfullplots
mv bias_${TAG}_test.png biases/bias_${TAG}.png
