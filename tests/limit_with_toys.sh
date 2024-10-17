#! /bin/bash
# Evaluate limits using toys
# Usage: ./limit_with_toys.sh <card>

Help() {
    echo "Evaluate limits with toys:"
    echo " Usage: limit_with_toys.sh <card> [options]"
    echo "--toys    (-t): N(toys) (default = 500)"
    echo "--gentoys (-g): N(gen) per segment (default = 100)"
    echo "--fitarg      : Additional fit arguments"
    echo "--genarg      : Additional generator arguments"
    echo "--name        : Name for output"
    echo "--tag         : Name tag for output"
    echo "--rrange  (-r): r range (default = 100)"
    echo "--seed    (-s): Base random seed (default = 90)"
    echo "--skipfits    : Skip fit loops, only create plots"
    echo "--grid        : Use a grid scan with MultiDimFit"
    echo "--dontclean   : Don't clean up temporary files"
}

CARD=""
NTOYS="500"
NGENPERTOY="100"
FITARG=""
GENARG=""
TAG=""
RRANGE="100"
SEED="90"
SKIPFITS=""
GRID=""
DONTCLEAN=""
NAME="test"

iarg=1
while [ ${iarg} -le $# ]
do
    eval "var=\${${iarg}}"
    if [[ "${var}" == "--help" ]] || [[ "${var}" == "-h" ]]
    then
        Help
        exit
    elif [[ "${var}" == "--rrange" ]] || [[ "${var}" == "-r" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        RRANGE=${var}
    elif [[ "${var}" == "--toys" ]] || [[ "${var}" == "-t" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        NTOYS=${var}
    elif [[ "${var}" == "--gentoys" ]] || [[ "${var}" == "-g" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        NGENPERTOY=${var}
    elif [[ "${var}" == "--name" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        NAME=${var}
    elif [[ "${var}" == "--tag" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        TAG=${var}
    elif [[ "${var}" == "--fitarg" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        FITARG=${var}
    elif [[ "${var}" == "--genarg" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        GENARG=${var}
    elif [[ "${var}" == "--seed" ]] || [[ "${var}" == "-s" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        SEED=${var}
    elif [[ "${var}" == "--skipfits" ]]
    then
        SKIPFITS="d"
    elif [[ "${var}" == "--grid" ]]
    then
        GRID="d"
    elif [[ "${var}" == "--dontclean" ]]
    then
        DONTCLEAN="d"
    elif [[ "${var}" == "--dryrun" ]]
    then
        DRYRUN="d"
    elif [[ "${CARD}" != "" ]]
    then
        echo "Arguments aren't configured correctly! (CARD = ${CARD}, var = ${var})"
        Help
        exit
    else
        CARD=${var}
    fi
    iarg=$((iarg + 1))
done

if [[ "${CARD}" == "" ]]; then
    Help
    exit
fi

if [ ! -f ${CARD} ]; then
    echo "No input card ${CARD} found"
    Help
    exit
fi

ARGUMENTS="--cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
ARGUMENTS="${ARGUMENTS} --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstant --X-rtd MINIMIZER_multiMin_maskConstraints"
ARGUMENTS="${ARGUMENTS} --rMin 0 --rMax ${RRANGE}"
ARGUMENTS=" --X-rtd MINIMIZER_multiMin_maskChannels=2"
ARGUMENTS="${ARGUMENTS} --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1"
if [[ "${CARD}" != *"total"* ]]; then
    ARGUMENTS="${ARGUMENTS} --cminRunAllDiscreteCombinations"
fi
ARGUMENTS="${ARGUMENTS} --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001"

combine -d ${CARD} -M HybridNew --rule CLs --cl=0.95 --LHCmode LHC-limits -T ${NTOYS} ${ARGUMENTS} ${FITARG}
