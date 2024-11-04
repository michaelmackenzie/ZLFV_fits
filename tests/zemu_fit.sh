#! /bin/bash
# Run a Z->emu fit
# Usage: ./zemu_fit.sh <card> [fit arg]

Help() {
    echo "Run a Z->emu fit"
    echo "Usage: ./zemu_fit.sh <card> [fit arg] [name] [r range]"
}

CARD=$1
FITARG=$2
NAME=$3
RRANGE=$4

if [[ "${CARD}" == "" ]]; then
    Help
    exit
fi

if [ ! -f ${CARD} ]; then
    echo "No input card ${CARD} found"
    Help
    exit
fi


if [[ "${RRANGE}" == "" ]]; then
    RMIN=-10
    RMAX=10
else
    RMIN=-${RRANGE}
    RMAX=${RRANGE}
fi

ARGUMENTS="--cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
ARGUMENTS="${ARGUMENTS} --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskConstraints" # --X-rtd MINIMIZER_multiMin_hideConstants
ARGUMENTS="${ARGUMENTS} --X-rtd MINIMIZER_multiMin_maskChannels=2"
ARGUMENTS="${ARGUMENTS} --rMin ${RMIN} --rMax ${RMAX}"

if [[ "${NAME}" == "" ]]; then
    NAME="test"
    if [[ "${CARD}" == *"bin1"* ]]; then
        NAME="bin1"
    elif [[ "${CARD}" == *"bin2"* ]]; then
        NAME="bin2"
    elif [[ "${CARD}" == *"bin3"* ]]; then
        NAME="bin3"
    elif [[ "${CARD}" == *"total"* ]]; then
        NAME="total"
    fi
fi
ARGUMENTS="-n .${NAME} ${ARGUMENTS}"

ARGUMENTS="${ARGUMENTS} --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1"
if [[ "${CARD}" != *"total"* ]]; then
    ARGUMENTS="${ARGUMENTS} --cminRunAllDiscreteCombinations"
fi
ARGUMENTS="${ARGUMENTS} --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001"

echo "Performing the fit"
echo ">>> combine -d ${CARD} -M FitDiagnostics ${ARGUMENTS} ${FITARG}"
combine -d ${CARD} -M FitDiagnostics ${ARGUMENTS} ${FITARG}

echo "Finished processing!"
