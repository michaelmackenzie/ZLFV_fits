#! /bin/bash

Help() {
    echo "Process combine impacts"
    echo "Usage: impacts.sh <card> [options]"
    echo "Options:"
    echo " --rrange    (-r ): POI range (default = 20)"
    echo " --obs       (-o ): Do observed"
    echo " --unblind        : Unblind the results"
    echo " --dontclean (-dc): Don't cleanup output combine files"
    echo " --approx    (-a ): Do approximate limits"
    echo " --plotarg        : Additional plotImpacts.py arguments (e.g. --blind, see plotImpacts.py -h)"
    echo " --fitarg         : Additional combineTool.py arguments (see combineTools.py -M Impacts -h)"
    echo " --summary        : Add summary info in the final plot"
    echo " --include   (-i ): Nuisance parameters to include (default is all)"
    echo " --exclude   (-e ): Exclude nuisance parameters (default is none)"
    echo " --fitonly        : Only process fitting steps, not impact PDF"
    echo " --skipinitial    : Skip initial fit"
    echo " --initialonly    : Only do initial fit"
    echo " --timeout        : Apply timeout to retry the initial fit"
    echo " --plotonly       : Skip fits, only make plots"
    echo " --parallel N     : Fit with N parallel processes"
    echo " --grid-sub       : Prepare for grid processing of impacts"
    echo " --grid-fin       : Finish grid processing of impacts"
    echo " --tag            : Tag for output results"
    echo " --dryrun         : Only print commands without processing"
    echo " --help      (-h ): Print this information"
}

CARD=""
DONTCLEAN=""
DOOBS=""
RRANGE="20"
APPROX=""
UNBLIND=""
COMMAND=""
PLOTARG=""
SUMMARY=""
EXCLUDE=""
INCLUDE=""
FITONLY=""
SKIPINITIAL=""
INITIALONLY=""
TIMEOUT=""
PLOTONLY=""
PARALLEL=""
GRIDSUB=""
GRIDFIN=""
TAG=""
VERBOSE=0
DRYRUN=""

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
    elif [[ "${var}" == "--exclude" ]] || [[ "${var}" == "-e" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        EXCLUDE=${var}
    elif [[ "${var}" == "--include" ]] || [[ "${var}" == "-i" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        INCLUDE=${var}
    elif [[ "${var}" == "--fitarg" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        COMMAND=${var}
    elif [[ "${var}" == "--plotarg" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        PLOTARG=${var}
    elif [[ "${var}" == "--approx" ]] || [[ "${var}" == "-a" ]]
    then
        APPROX="d"
    elif [[ "${var}" == "--obs" ]] || [[ "${var}" == "-o" ]]
    then
        DOOBS="d"
    elif [[ "${var}" == "--unblind" ]]
    then
        UNBLIND="d"
    elif [[ "${var}" == "--fitonly" ]]
    then
        FITONLY="d"
    elif [[ "${var}" == "--summary" ]]
    then
        SUMMARY="d"
    elif [[ "${var}" == "--skipinitial" ]]
    then
        SKIPINITIAL="d"
    elif [[ "${var}" == "--initialonly" ]]
    then
        INITIALONLY="d"
    elif [[ "${var}" == "--timeout" ]]
    then
        TIMEOUT="d"
    elif [[ "${var}" == "--plotonly" ]]
    then
        PLOTONLY="d"
    elif [[ "${var}" == "--parallel" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        PARALLEL=${var}
    elif [[ "${var}" == "--grid-sub" ]]
    then
        GRIDSUB="d"
    elif [[ "${var}" == "--grid-fin" ]]
    then
        GRIDFIN="d"
    elif [[ "${var}" == "--grid" ]]
    then
        echo "--grid isn't defined, use --grid-sub for grid submission and --grid-fin for finishing a grid run"
        Help
        exit
    elif [[ "${var}" == "--tag" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        TAG=${var}
    elif [[ "${var}" == "--dontclean" ]] || [[ "${var}" == "-dc" ]]
    then
        DONTCLEAN="d"
    elif [[ "${var}" == "--verbose" ]] || [[ "${var}" == "-v" ]]
    then
        #Test if a verbosity level is given
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        if [[ $var =~ ^-?[0-9]+$ ]]
        then
            VERBOSE=${var}
        else
            iarg=$((iarg - 1))
            eval "var=\${${iarg}}"
            VERBOSE=1
        fi
    elif [[ "${var}" == "--dryrun" ]]
    then
        DRYRUN="d"
    elif [[ "${CARD}" != "" ]]
    then
        echo "Arguments aren't configured correctly! card = ${CARD}, arg = ${var} "
        echo $@
        Help
        exit
    else
        CARD=${var}
    fi
    iarg=$((iarg + 1))
done

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]
then
    Help
    exit
fi

if [[ "${CARD}" == "" ]]
then
    echo "No combine card given!"
    Help
    exit
fi

if [ ! -f ${CARD} ]
then
    echo "Combine card ${CARD} not found"
    Help
    exit
fi

if [[ "${DOOBS}" == "" ]]
then
    DOOBS="-t -1"
else
    DOOBS=""
fi

if [[ "${RRANGE}" == "" ]]
then
    RRANGE=20
fi

if [[ "${GRIDFIN}" != "" ]]; then
    DONTCLEAN="d"
    SKIPINITIAL="d"
fi

echo "Performing impacts on combine card ${CARD}"
echo "Parameters: DONTCLEAN=${DONTCLEAN}; DOOBS=${DOOBS}; RRANGE=${RRANGE}"

WORKSPACE=`echo ${CARD} | sed 's/.txt/_workspace.root/'`
JSON=`echo ${CARD} | sed 's/.txt/.json/' | sed 's/_workspace.root/.json/' | sed 's/.root/.json/' | sed 's/datacard_/impacts_/'`
FITNAME=`echo ${CARD} | sed 's/.txt//' | sed 's/_workspace.root//' | sed 's/.root//' | sed 's/datacard_zprime_//'`
# Task name for grid processing
TASK=`echo ${CARD} | sed 's/.txt//' | sed 's/_workspace.root//' | sed 's/.root//'`
if [[ "${APPROX}" != "" ]]; then
    JSON=`echo ${JSON} | sed "s/.json/_approx.json/"`
fi
if [[ "${TAG}" != "" ]]; then
    JSON=`echo ${JSON} | sed "s/.json/_${TAG}.json/"`
    FITNAME="${FITNAME}_${TAG}"
    TASK="${TASK}_${TAG}"
fi

echo "FITNAME=${FITNAME}"


if [[ "${PLOTONLY}" == "" ]]; then
    ARGS=""
    if [[ "${APPROX}" != "" ]]; then
        echo "-- Performing approximate impacts"
        ARGS="--approx robust"
    fi
    if [[ "${COMMAND}" != "" ]]; then
        echo "Adding argument ${COMMAND} to the evaluation"
        ARGS="${ARGS} ${COMMAND}"
    fi

    if [[ "${CARD}" == *".txt" ]] && [[ "${GRIDFIN}" == "" ]]; then
        COMMAND="text2workspace.py ${CARD} -o ${WORKSPACE}"
        echo ${COMMAND}
        if [[ "${DRYRUN}" == "" ]]; then
            ulimit -s unlimited
            ${COMMAND}
        fi
    fi

    FITADDITIONAL=""
    if [[ "${EXCLUDE}" != "" ]]; then
        FITADDITIONAL="${FITADDITONAL} --exclude ${EXCLUDE}"
    fi
    if [[ "${INCLUDE}" != "" ]]; then
        FITADDITIONAL="${FITADDITONAL} --named ${INCLUDE}"
    fi

    # Add additional fit arguments to help convergence
    FITADDITIONAL="${FITADDITIONAL} --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
    FITADDITIONAL="${FITADDITIONAL} --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskConstraints"
    FITADDITIONAL="${FITADDITIONAL} --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --X-rtd MINIMIZER_multiMin_maskChannels=2"
    FITADDITIONAL="${FITADDITIONAL} --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001"

    # Perform initial fits
    if [[ "${APPROX}" == "" ]] && [[ "${SKIPINITIAL}" == "" ]] && [[ "${GRIDFIN}" == "" ]]; then
        COMMAND="combineTool.py -M Impacts -d ${WORKSPACE} -m 0 --rMin -${RRANGE} --rMax ${RRANGE} --robustFit 1 --doInitialFit ${DOOBS} ${FITADDITIONAL} -n ${FITNAME}"
        if [[ "${TIMEOUT}" != "" ]]; then
            COMMAND="timeout 300 ${COMMAND}"
        fi
        echo ${COMMAND}
        if [[ "${DRYRUN}" == "" ]]; then
            ${COMMAND}
            if [ $? -ne 0 ]; then
                echo "Attempting initial fit again due to failure/timeout!"
                ${COMMAND}
            fi
        fi
    fi

    if [[ "${INITIALONLY}" != "" ]]; then
        echo "Finished initial processing"
        exit
    fi

    # Perform impact fits
    COMMAND="combineTool.py -M Impacts -d ${WORKSPACE} -m 0 --rMin -${RRANGE} --rMax ${RRANGE} --robustFit 1 --doFits ${DOOBS} ${ARGS} ${FITADDITIONAL} -n ${FITNAME}"
    if [[ "${PARALLEL}" != "" ]]; then
        COMMAND="${COMMAND} --parallel ${PARALLEL}"
    fi
    if [[ "${GRIDSUB}" != "" ]]; then
        echo "Using grid task name ${TASK}"
        RESULTDIR=grid/${TASK}
        [ ! -f ${RESULTDIR} ] && mkdir -p ${RESULTDIR}
        cd ${RESULTDIR}
        ln -s ../../${WORKSPACE} ${WORKSPACE}
        mv ../../higgsCombine_initialFit_${FITNAME}.MultiDimFit.mH0.root ./
        echo ">>> ${COMMAND} --job-mode condor --task-name ${TASK} --sub-opts +DesiredOS = \"SL7\""
        if [[ "${DRYRUN}" == "" ]]; then
            ${COMMAND} --job-mode condor --task-name ${TASK} --sub-opts '+DesiredOS = "SL7"'
        fi
        TASKDIR=~/private/grid/combine/${TASK}
        [ -d ${TASKDIR} ] && rm ${TASKDIR}/*
        [ ! -d ${TASKDIR} ] && mkdir -p ${TASKDIR}
        COMMAND="mv condor_${TASK}.s* ${TASKDIR}/"
        echo ">>> ${COMMAND}"
        if [[ "${DRYRUN}" == "" ]]; then
            ${COMMAND}
        fi
        echo "Exiting with grid setup complete, use grid directory ${TASKDIR}"
        exit
    fi

    # Local processing

    if [[ "${GRIDFIN}" == "" ]]; then
        if [[ "${TIMEOUT}" != "" ]]; then
            COMMAND="timeout 500 ${COMMAND}"
        fi
        echo ${COMMAND}
        if [[ "${DRYRUN}" == "" ]]; then
            ${COMMAND}
            if [ $? -ne 0 ]; then
                echo "Attempting initial fit again due to failure/timeout!"
                ${COMMAND}
            fi
        fi
    else
        RESULTDIR=grid/${TASK}
        cd ${RESULTDIR}
    fi

    # Compile impact fits

    COMMAND="combineTool.py -M Impacts -d ${WORKSPACE} -m 0 --rMin -${RRANGE} --rMax ${RRANGE} --robustFit 1 --output ${JSON} ${DOOBS} ${ARGS} ${FITADDITIONAL} -n ${FITNAME}"
    echo ">>> ${COMMAND}"
    if [[ "${DRYRUN}" == "" ]]; then
        ${COMMAND}
    fi
fi

if [[ "${FITONLY}" ]]; then
    exit
fi

if [[ "${PLOTONLY}" != "" ]] && [[ "${GRIDFIN}" != "" ]]; then
    RESULTDIR=grid/${TASK}
    cd ${RESULTDIR}
fi


#Blind the result if using the observed impacts but not yet unblinding the measurement
# DOOBS = -t -1 == Asimov
if [[ "${DOOBS}" == "" ]] && [[ "${UNBLIND}" == "" ]]; then
    PLOTARG="${PLOTARG} --blind"
fi
if [[ "${SUMMARY}" != "" ]]; then
    PLOTARG="${PLOTARG} --summary"
fi

COMMAND="plotImpacts.py -i ${JSON} ${PLOTARG} -o `echo ${JSON} | sed 's/.json//'`"
echo ">>> ${COMMAND}"
if [[ "${DRYRUN}" == "" ]]; then
    ${COMMAND}
fi

if [[ "${DONTCLEAN}" == "" ]]
then
    COMMAND="rm higgsCombine_*_${FITNAME}*.MultiDimFit.mH0.root"
    echo ">>> ${COMMAND}"
    if [[ "${DRYRUN}" == "" ]]; then
        ${COMMAND}
    fi
fi
