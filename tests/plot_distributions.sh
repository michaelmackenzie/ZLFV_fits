#! /bin/bash
# Create pre- and post-fit distributions

Help() {
    echo "Plot pre- and post-fit distributions"
    echo "Usage: plot_distributions.sh <card> [options]"
    echo "Options:"
    echo " --rrange    (-r ): POI range (default = 30)"
    echo " --unblind   (-u ): Unblind post-fit distributions"
    echo " --fitarg         : Additional fit arguments"
    echo " --ignoresys      : Don't save shape uncertainties in output (for debugging)"
    echo " --skipworkspace  : Don't recreate workspace for input data card"
    echo " --plotonly       : Only make plots"
    echo " --plotparams     : Make plots of all fit parameters"
    echo " --zemu           : Assume Z->e+mu processing"
    echo " --tag            : Tag for output results"
    echo " --help      (-h ): Print this information"
}

CARD=""
RRANGE="30"
UNBLIND="false"
FITARG=""
IGNORESYS=""
SKIPWORKSPACE=""
PLOTONLY=""
PLOTPARAMS=""
ZEMU=""
TAG=""

iarg=1
while [ ${iarg} -le $# ]
do
    eval "var=\${${iarg}}"
    if [[ "${var}" == "--help" ]] || [[ "${var}" == "-h" ]]; then
        Help
        exit
    elif [[ "${var}" == "--rrange" ]] || [[ "${var}" == "-r" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        RRANGE=${var}
    elif [[ "${var}" == "--unblind" ]] || [[ "${var}" == "-u" ]]; then
        UNBLIND="true"
    elif [[ "${var}" == "--fitarg" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        FITARG=${var}
    elif [[ "${var}" == "--tag" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        TAG=${var}
    elif [[ "${var}" == "--ignoresys" ]]; then
        IGNORESYS="d"
    elif [[ "${var}" == "--skipworkspace" ]]; then
        SKIPWORKSPACE="d"
    elif [[ "${var}" == "--plotonly" ]]; then
        PLOTONLY="d"
    elif [[ "${var}" == "--plotparams" ]]; then
        PLOTPARAMS="d"
    elif [[ "${var}" == "--zemu" ]]; then
        ZEMU="d"
    elif [[ "${CARD}" != "" ]]; then
        echo "Arguments aren't configured correctly! (at ${var}, CARD=${CARD})"
        Help
        exit
    else
        CARD=${var}
    fi
    iarg=$((iarg + 1))
done

if [[ "${CARD}" == "" ]]; then
    echo "No card given!"
    exit
fi

FIGUREPATH="figures"
if [[ "${TAG}" != "" ]]; then
    FIGUREPATH="${FIGUREPATH}/${TAG}"
    TAG="_${TAG}"
fi
echo "Using card ${CARD} with output figure path ${FIGUREPATH}"

# Create a workspace if needed

WORKSPACE=`echo ${CARD} | sed 's/.txt/_workspace.root/'`
if [[ "${SKIPWORKSPACE}" == "" ]] && [[ "${PLOTONLY}" == "" ]] && [[ "${CARD}" == *".txt" ]]; then
    echo "Creating workspace ${WORKSPACE}"
    text2workspace.py ${CARD} --channel-masks -o ${WORKSPACE}
    echo "Created workspace ${WORKSPACE}"
fi

FITARG="${FITARG} --saveShapes --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1 --rMin -${RRANGE} --rMax ${RRANGE}"
if [[ "${ZEMU}" != "" ]]; then
    FITARG="${FITARG} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams" #--X-rtd MINIMIZER_multiMin_hideConstants
    FITARG="${FITARG} --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
    FITARG="${FITARG} --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001"
    if [[ "${WORKSPACE}" != *"total"* ]] && [[ "${WORKSPACE}" != *"zprime"* ]]; then
        FITARG="${FITARG} --cminRunAllDiscreteCombinations"
    fi
fi
if [[ "${IGNORESYS}" == "" ]]; then
    FITARG="${FITARG} --saveWithUncertainties"
fi

if [[ "${PLOTONLY}" == "" ]]; then
    echo "Performing the initial fit"
    echo ">>> combine -M FitDiagnostics -d ${WORKSPACE} ${FITARG} -n .plot_distributions${TAG}"
    combine -M FitDiagnostics -d ${WORKSPACE} ${FITARG} -n .plot_distributions${TAG}
fi

FILE="fitDiagnostics.plot_distributions${TAG}.root"
if [ ! -f ${FILE} ]; then
    echo "${FILE} not found!"
    exit
fi

if [[ "${ZEMU}" == "" ]]; then
    root.exe -q -b "${CMSSW_BASE}/src/ZLFV_fits/tools/plot_fit.C(\"${FILE}\", \"${FIGUREPATH}\", ${UNBLIND})"
else
    root.exe -q -b "${CMSSW_BASE}/src/ZLFV_fits/tools/plot_bemu.C(\"${FILE}\", \"${FIGUREPATH}\", ${UNBLIND})"
fi

if [[ "${PLOTPARAMS}" ]]; then
    root.exe -q -b "${CMSSW_BASE}/src/ZLFV_fits/tools/plot_combine_fit_params.C(\"${FILE}\", \"${FIGUREPATH}_params\")"
fi
echo "Finished running distribution plotting"
