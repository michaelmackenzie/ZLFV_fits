#! /bin/bash
# Test the look elsewhere effect (LEE) for a resonance scan by fitting toys

Help() {
    echo "Evaluate the look elsewhere effect for a resonance scan:"
    echo " Usage: bemu_lee.sh <scan base name> [options]"
    echo "--toys    (-t): N(toys) (default = 1000)"
    echo "--name        : Name for output"
    echo "--tag         : Name tag for output"
    echo "--makebkgpdf  : Create the full data spectrum PDFs for toy generation"
    echo "--skiptoygen  : Skip creation of the toys/datacards"
    echo "--skipfits    : Skip fit loops, only create plots"
    echo "--skipplots   : Skip plot creation"
    echo "--dontclean   : Don't clean up temporary files"
}

BASENAME=""
NTOYS=100
TAG=""
MAKEBKGPDF=""
SKIPTOYGEN=""
SKIPFITS=""
SKIPPLOTS=""
DONTCLEAN=""
DRYRUN=""

iarg=1
while [ ${iarg} -le $# ]
do
    eval "var=\${${iarg}}"
    if [[ "${var}" == "--help" ]] || [[ "${var}" == "-h" ]]
    then
        Help
        exit
    elif [[ "${var}" == "--toys" ]] || [[ "${var}" == "-t" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        NTOYS=${var}
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
    elif [[ "${var}" == "--makebkgpdf" ]]
    then
        MAKEBKGPDF="d"
    elif [[ "${var}" == "--skiptoygen" ]]
    then
        SKIPTOYGEN="d"
    elif [[ "${var}" == "--skipfits" ]]
    then
        SKIPFITS="d"
    elif [[ "${var}" == "--skipplots" ]]
    then
        SKIPPLOTS="d"
    elif [[ "${var}" == "--dontclean" ]]
    then
        DONTCLEAN="d"
    elif [[ "${var}" == "--dryrun" ]]
    then
        DRYRUN="d"
    elif [[ "${BASENAME}" != "" ]]
    then
        echo "Arguments aren't configured correctly! (BASENAME = ${BASENAME}, var = ${var})"
        Help
        exit
    else
        BASENAME=${var}
    fi
    iarg=$((iarg + 1))
done

if [[ "${BASENAME}" == "" ]]; then
    Help
    exit
fi

# Re-create the full dataset PDF
if [[ "${MAKEBKGPDF}" != "" ]]; then
    python ScanMuE_fit_wrapper_v2.py -o bdt_0d3_0d7_LEE --full-mass --scan-min 300 --scan-max 300.1 --scan-step 1 --xgb-min 0.30 --xgb-max 0.70 --param-name bin1 --component bkg
    python ScanMuE_fit_wrapper_v2.py -o bdt_0d7_1d0_LEE --full-mass --scan-min 300 --scan-max 300.1 --scan-step 1 --xgb-min 0.70 --xgb-max 1.01 --param-name bin2 --component bkg
fi

# Generate toy datasets to use
if [[ "${SKIPTOYGEN}" == "" ]]; then

    # Create each toy dataset
    COUNT=0 # For processing data card creation in parallel
    [ ! -d log ] && mkdir log
    for ((ITOY=1; ITOY<=NTOYS;ITOY++)); do
        echo "Creating toy ${ITOY}..."
        ((COUNT++))
        python create_toy.py --fit-file WorkspaceScanBKG/workspace_scanbkg_v2_bdt_0d3_0d7_LEE_mp0.root -o lee_toy_0d3_0d7 --toy ${ITOY} --param bin1 --seed ${ITOY}
        python create_toy.py --fit-file WorkspaceScanBKG/workspace_scanbkg_v2_bdt_0d7_1d0_LEE_mp0.root -o lee_toy_0d7_1d0 --toy ${ITOY} --param bin2 --seed ${ITOY}
        echo "Creating the corresponding datacards..."
        LOGFILE="log/lee_card_creation_${ITOY}.log"
        if [ -f ${LOGFILE} ]; then
            rm ${LOGFILE}
        fi
        ./clone_cards_for_toy.sh datacards/${BASENAME}/ datacards/${BASENAME}_toy_${ITOY}/ WorkspaceScanTOY/toy_file_lee_toy_0d3_0d7_${ITOY}.root WorkspaceScanTOY/toy_file_lee_toy_0d7_1d0_${ITOY}.root >| ${LOGFILE} 2>&1 &
        if [ ${COUNT} -ge 8 ]; then
            COUNT=0
            wait
        fi
    done
    wait
fi

# Fit each toy
if [[ "${SKIPFITS}" == "" ]]; then
    for ((ITOY=1; ITOY<=NTOYS;ITOY++)); do
        echo "Fitting toy ${ITOY}..."
        python perform_scan.py -o ${BASENAME}_toy_${ITOY} --unblind --smooth-expected
    done
fi

# Create a LEE map
python tools/toys_to_lee_map.py -o ${BASENAME} --ntoys ${NTOYS} --verbose 1
