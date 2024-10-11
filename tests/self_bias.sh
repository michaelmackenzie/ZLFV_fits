#! /bin/bash
# Test the bias when generating with a given PDF and fitting with the envelope

CARD=$1
BIN=$2
NPDFS=$3
TAG=$4
TOYS=$5

if [[ "${CARD}" == "" ]]; then
    echo "No input card provided!"
    exit
fi
if [ ! -f ${CARD} ]; then
    echo "Input card ${CARD} not found!"
    exit
fi
if [[ "${NPDFS}" == "" ]]; then
    NPDFS=6
fi
if [[ "${BIN}" == "" ]]; then
    BIN="bin1"
fi
if [[ "${TOYS}" == "" ]]; then
    TOYS="500"
fi

for (( IPDF=0; IPDF<${NPDFS}; IPDF++ )); do
    ARGS="--card_1 ${CARD} --card_2 ${CARD} --name ${TAG}_${IPDF} -g 500 -t ${TOYS} -r 50 --skipfullplots"
    ${CMSSW_BASE}/src/ZLFV_fits/tests/bemu_gen_fit_test.sh ${ARGS} --genarg "--setParameters pdfindex_${BIN}=${IPDF}"
    mv bias_${TAG}_${IPDF}_test.png biases/bias_${TAG}_${IPDF}.png
done
