#! /bin/bash

CARD=$1
VAR=$2
VAL=$3
RPOINT=$4


ARGUMENTS="--cminDefaultMinimizerStrategy 0 --saveNLL -m 125 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
ARGUMENTS="${ARGUMENTS} --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstant --X-rtd MINIMIZER_multiMin_maskConstraints"
ARGUMENTS="${ARGUMENTS} --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --X-rtd MINIMIZER_multiMin_maskChannels=2"
ARGUMENTS="${ARGUMENTS} --setParameters ${VAR}=${VAL},r=${RPOINT} --freezeParameters r,${VAR}"

combine -M MultiDimFit -d ${CARD} ${ARGUMENTS} -v 3 --rMin -5 --rMax 5
