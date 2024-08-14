#!/bin/bash

CODE=$1
shift;

PROCESS=$1
shift;

NAME=$1
shift;


#FLAGS="-t -1 --expectSignal 2"
#FLAGS="-t -1 --expectSignal 0"
FLAGS=""

#FROZENPARAM="--freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3 --setParameters pdfindex_bin1=0,pdfindex_bin2=0,pdfindex_bin3=0"
FROZENPARAM=""

TRANSFER=1
CLEANCOMB=0
ORGDIR=$(pwd)



COMBpath=/afs/cern.ch/work/g/gkaratha/private/Analysis/DispJets/Analyzer/Limit/TutorialCombine/CMSSW_11_3_4/src/combinetutorial-2023-parametric/

cp $CODE $COMBpath
cp workspace*.root $COMBpath


cd $COMBpath
cmsenv

echo ">>>>>> process ${PROCESS} with flags: $FLAGS AND fixed parameters: $FROZENPARAM <<<<<<<<<"


text2workspace.py $CODE -m 91 -o param_datacard_${NAME}.root


if  [[ "${PROCESS}" == "asimovlim" ]]
then
   combine -M AsymptoticLimits param_datacard_${NAME}.root  -m 91 --saveWorkspace -n .bestfit_${NAME} -t -1 $FROZENPARAM  $FLAGS
fi

if  [[ "${PROCESS}" == "limits" ]]
then
   combine -M AsymptoticLimits param_datacard_${NAME}.root  -m 91 --saveWorkspace -n .bestfit_${NAME} $FROZENPARAM  $FLAGS
fi

if  [[ "${PROCESS}" == "fit" ]]
then
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 --saveWorkspace -n .bestfit_${NAME} $FROZENPARAM  $FLAGS --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2
   echo "higgsCombine.bestfit_${NAME}.MultiDimFit.mH91.root has been created in Roofit folder - No need to move"
   echo "for plot run:"
   echo  "python multdimfit_plotter.py -f higgsCombine.bestfit_${NAME}.MultiDimFit.mH91.root -o ${NAME} --bin <bin1> --mllvar <lepm_11> (-v)" 
fi

if  [[ "${PROCESS}" == "conf" ]]
then
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 --saveWorkspace -n .singles_${NAME} --algo singles $FROZENPARAM  $FLAGS
fi


if  [[ "${PROCESS}" == "scan" ]]
then
#   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .scan_r0_${NAME} --algo grid --points 10 --setParameterRanges r=-5,5 -t -1 --expectSignal 0
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .scan_r0_${NAME} --algo grid --points 20 --setParameterRanges r=-5.0,5.0 -t -1 --expectSignal 0
   plot1DScan.py higgsCombine.scan_r0_${NAME}.MultiDimFit.mH91.root -o ${NAME}_scan_r0
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .scan_r1_${NAME} --algo grid --points 20 --setParameterRanges r=-5.0,5.0 -t -1 --expectSignal 1
   plot1DScan.py higgsCombine.scan_r1_${NAME}.MultiDimFit.mH91.root -o ${NAME}_scan_r1
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .scan_r2_${NAME} --algo grid --points 20 --setParameterRanges r=-5.0,5.0 -t -1 --expectSignal 2
   plot1DScan.py higgsCombine.scan_r2_${NAME}.MultiDimFit.mH91.root -o ${NAME}_scan_r2
   mv ${NAME}_scan_r*.png $ORGDIR
   echo "${NAME}_scan has been created"
fi


if  [[ "${PROCESS}" == "gof" ]]
then
   echo "runs ONLY on binned fits - needs real data"
   text2workspace.py $CODE -m 91 -o param_datacard_${NAME}.root --channel-masks
   combine -M GoodnessOfFit param_datacard_${NAME}.root --algo saturated -m 91.0 -n .goodnessOfFit_data_${NAME} $FROZENPARAM
   combine -M GoodnessOfFit param_datacard_${NAME}.root --algo saturated -m 91.0 -n .goodnessOfFit_toys_${NAME} -t 300 $FROZENPARAM
   combineTool.py -M CollectGoodnessOfFit --input higgsCombine.goodnessOfFit_data_${NAME}.GoodnessOfFit.mH91.root higgsCombine.goodnessOfFit_toys_${NAME}.GoodnessOfFit.mH91.123456.root -m 91.0 -o gof_${NAME}.json
   plotGof.py gof_${NAME}.json --statistic saturated --mass 91.0 -o ${NAME}_gof
   mv ${NAME}_gof.png $ORGDIR
   mv ${NAME}_gof.pdf $ORGDIR
   echo "${NAME}_gof.png /pdf has been created"
fi


if  [[ "${PROCESS}" == "unc" ]]
then
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .scan.with_syst_${NAME} --algo grid --points 20 --setParameterRanges r=-5.0,5.0 -t -1 --expectSignal 0
   combine -M MultiDimFit param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .bestfit.with_syst_${NAME} --saveWorkspace -t -1 --expectSignal 0 --setParameterRanges r=-5.0,5.0
 if [[ "${FROZENPARAM}" == "" ]]
    then
      FROZENPARAM="--freezeParameters "
    fi
   combine -M MultiDimFit  higgsCombine.bestfit.with_syst_${NAME}.MultiDimFit.mH91.root -m 91 ${FROZENPARAM},allConstrainedNuisances -n .scan.statonly_${NAME} --algo grid --points 20 --snapshotName MultiDimFit -t -1 --expectSignal 0 --setParameterRanges r=-5.0,5.0
   plot1DScan.py higgsCombine.scan.with_syst_${NAME}.MultiDimFit.mH91.root --main-label Total --main-color 1 --others higgsCombine.scan.statonly_${NAME}.MultiDimFit.mH91.root:"Stat-only":2 -o ${NAME}_unc_scan --breakdown exp,stat
   mv ${NAME}_unc_scan.png $ORGDIR
   echo "${NAME}_unc_scan has been created"
fi

if  [[ "${PROCESS}" == "impacts" ]]
then
   EXTRA_OPT=" --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --robustFit 1 "
   #--cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --robustFit 1
   #"--freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3 --setParameters pdfindex_bin1=0,pdfindex_bin2=0,pdfindex_bin3=0"
   #--cminRunAllDiscreteCombinations
   # --setParameterRanges r=-1,1
   #--X-rtd MINIMIZER_multiMin_maskChannels= 1 or 2
   # --setParameterRanges shapeBkg_bkg_bin1__norm=24000,25500:shapeBkg_bkg_bin2__norm=22000,30000
   # --rMin -10 --rMax 10
   combineTool.py -M Impacts -d param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .impacts_${NAME}  --doInitialFit $EXTRA_OPT 
   combineTool.py -M Impacts -d param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .impacts_${NAME} --doFits $EXTRA_OPT --parallel 10
   combineTool.py -M Impacts -d param_datacard_${NAME}.root -m 91 $FROZENPARAM -n .impacts_${NAME} -o impacts_${NAME}.json $EXTRA_OPT
   plotImpacts.py -i impacts_${NAME}.json -o ${NAME}_impacts --blind --summary
   mv ${NAME}_impacts.pdf $ORGDIR
   echo "${NAME}_impacts.pdf has been created"
fi


if [[ "${PROCESS}" == "closure-envelope" ]]
then
    
    GEN_CODE=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;

    cp $ORGDIR/$GEN_CODE .
    echo "Run closure test with template $GEN_CODE, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"

    text2workspace.py $GEN_CODE -m 91 -o param_template_${NAME}.root    
    combine -M GenerateOnly param_template_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .generate_template_${NAME}  --expectSignal $R_FOR_BIAS --saveToys 


    combine -M MultiDimFit param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .closure_envelope_${NAME}  --expectSignal $R_FOR_BIAS --toysFile higgsCombine.generate_template_${NAME}.GenerateOnly.mH91.123456.root --algo singles  --rMin -100 --rMax 100 --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace --saveToys  #--cminRunAllDiscreteCombinations
    python3 plot_bias_pull.py -i higgsCombine.closure_envelope_${NAME}.MultiDimFit.mH91.123456.root -o ${NAME}_closure_envelope --n-toys $TOYS_FOR_BIAS --r-truth $R_FOR_BIAS --gen-fnc mc --fit-fnc all
   
    mv ${NAME}_closure_envelope_pull.png $ORGDIR
    echo "${NAME}_closure_envelope_pull.png has been created"
    mv ${NAME}_closure_envelope_pull_withfit.png $ORGDIR
    echo "${NAME}_closure_envelope_pull_withfit.png has been created"
   

fi


if [[ "${PROCESS}" == "closure-pdf" ]]
then
    
    GEN_CODE=$1
    shift;
    FIT_PDF_IDX=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;

    INDEX_FIT=""

    if [[ "${BRANCH}" == "" ]]
    then
       INDEX_FIT="--setParameters pdfindex_bin1=$FIT_PDF_IDX,pdfindex_bin2=$FIT_PDF_IDX,pdfindex_bin3=$FIT_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3 --saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else
       INDEX_FIT="--setParameters $BRANCH=$FIT_PDF_IDX --freezeParameters $BRANCH --saveSpecifiedIndex $BRANCH"
    fi


    cp $ORGDIR/$GEN_CODE .
    echo "Run closure test with template $GEN_CODE, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"

    text2workspace.py $GEN_CODE -m 91 -o param_template_${NAME}.root    
    combine -M GenerateOnly param_template_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .generate_template_${NAME} --expectSignal $R_FOR_BIAS --saveToys 

    
    combine -M MultiDimFit param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .closure_pdf_${FIT_PDF_IDX}_${NAME}  --expectSignal $R_FOR_BIAS --toysFile higgsCombine.generate_template_${NAME}.GenerateOnly.mH91.123456.root --algo singles  --rMin -100 --rMax 100 --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace --saveToys ${INDEX_FIT} 
    python3 plot_bias_pull.py -i higgsCombine.closure_pdf_${FIT_PDF_IDX}_${NAME}.MultiDimFit.mH91.123456.root -o ${NAME}_closure_pdf_${FIT_PDF_IDX} --n-toys $TOYS_FOR_BIAS --r-truth $R_FOR_BIAS --gen-fnc mc --fit-fnc ${FIT_PDF_IDX}
   
    mv ${NAME}_closure_pdf_${FIT_PDF_IDX}_pull.png $ORGDIR
    echo "${NAME}_closure_pdf_${FIT_PDF_IDX}_pull.png has been created"
    mv ${NAME}_closure_pdf_${FIT_PDF_IDX}_pull_withfit.png $ORGDIR
    echo "${NAME}_closure_pdf_${FIT_PDF_IDX}_pull_withfit.png has been created"

fi





if [[ "${PROCESS}" == "bias-pdf" ]]
then
    
    GEN_PDF_IDX=$1
    shift;
    FIT_PDF_IDX=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;

    INDEX_TOY=""
    INDEX_FIT=""

    if [[ "${BRANCH}" == "" ]]
    then
      INDEX_TOY="--setParameters pdfindex_bin1=$GEN_PDF_IDX,pdfindex_bin2=$GEN_PDF_IDX,pdfindex_bin3=$GEN_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
      INDEX_FIT="--setParameters pdfindex_bin1=$FIT_PDF_IDX,pdfindex_bin2=$FIT_PDF_IDX,pdfindex_bin3=$FIT_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3 --saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else 
       INDEX_TOY="--setParameters $BRANCH=$GEN_PDF_IDX --freezeParameters $BRANCH"
       INDEX_FIT="--setParameters $BRANCH=$FIT_PDF_IDX --freezeParameters $BRANCH --saveSpecifiedIndex $BRANCH"
    fi 

    echo "Run bias of idividual padfs on the envelope, pdf index for toys = $GEN_PDF_IDX, pdf index for fit = $FIT_PDF_IDX, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"
    echo $INDEX_TOY
    echo $INDEX_FIT


    combine -M GenerateOnly param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .generate_truth_index_${GEN_PDF_IDX}_${NAME} --expectSignal $R_FOR_BIAS --saveToys ${INDEX_TOY}
 
    combine -M MultiDimFit param_datacard_${NAME}.root --algo singles -m 91 -t $TOYS_FOR_BIAS -n .bias_truth_index_${GEN_PDF_IDX}_fit_index_${FIT_PDF_IDX}_${NAME} --toysFile higgsCombine.generate_truth_index_${GEN_PDF_IDX}_${NAME}.GenerateOnly.mH91.123456.root --rMin -100 --rMax 100 --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace --saveToys ${INDEX_FIT}

#--trackParameters cheb3_bkgPDF_bin1_0,cheb3_bkgPDF_bin1_1,cheb3_bkgPDF_bin1_2,cheb3_bkgPDF_bin1_0,cheb3_bkgPDF_bin1_1,sumexp2_bkgPDF_bin1_x0,sumexp2_bkgPDF_bin1_x1,sumexp2_bkgPDF_bin1_c0,sumexp2_bkgPDF_bin1_c1,pdfindex_bin1,shapeBkg_bkg_bin1__norm,shapeSig_sgn_bin1__norm

    python3 plot_bias_pull.py -i higgsCombine.bias_truth_index_${GEN_PDF_IDX}_fit_index_${FIT_PDF_IDX}_${NAME}.MultiDimFit.mH91.123456.root -o ${NAME}_truth_${GEN_PDF_IDX}_fit_${FIT_PDF_IDX} --n-toys $TOYS_FOR_BIAS --r-truth $R_FOR_BIAS --gen-fnc PDF_${GEN_PDF_IDX} --fit-fnc PDF_${FIT_PDF_IDX}
   
    mv ${NAME}_truth_${GEN_PDF_IDX}_fit_${FIT_PDF_IDX}_pull.png $ORGDIR
    echo "${NAME}_truth_${GEN_PDF_IDX}_fit_${FIT_PDF_IDX}_pull.png has been created"
    mv ${NAME}_truth_${GEN_PDF_IDX}_fit_${FIT_PDF_IDX}_pull_withfit.png $ORGDIR
    echo "${NAME}_truth_${GEN_PDF_IDX}_fit_${FIT_PDF_IDX}_pull_withfit.png has been created"
    

fi


if [[ "${PROCESS}" == "bias-envelope" ]]
then
    
    GEN_PDF_IDX=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;
    
    INDEX_TOY=""
    INDEX_FIT=""
    if [[ "${BRANCH}" == "" ]]
    then
      INDEX_TOY="--setParameters pdfindex_bin1=$GEN_PDF_IDX,pdfindex_bin2=$GEN_PDF_IDX,pdfindex_bin3=$GEN_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
      INDEX_FIT="--saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else 
       INDEX_TOY="--setParameters $BRANCH=$GEN_PDF_IDX --freezeParameters $BRANCH"
       INDEX_FIT="--saveSpecifiedIndex $BRANCH --cminRunAllDiscreteCombinations"
    fi 

    echo "Run bias for the envelope (gen pdf from envelope), pdf index for toys = $GEN_PDF_IDX, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"
    echo ${INDEX_TOY}


    combine -M GenerateOnly param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .generate_truth_index_${GEN_PDF_IDX}_${NAME}  --expectSignal $R_FOR_BIAS --saveToys  ${INDEX_TOY}
 
    combine -M MultiDimFit param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .bias_truth_index_${GEN_PDF_IDX}_fit_envelop_${NAME} --toysFile higgsCombine.generate_truth_index_${GEN_PDF_IDX}_${NAME}.GenerateOnly.mH91.123456.root --algo singles  --rMin -100 --rMax 100 --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace ${INDEX_FIT} # --cminRunAllDiscreteCombinations
#--trackParameters cheb3_bkgPDF_bin2_0,cheb3_bkgPDF_bin2_1,cheb3_bkgPDF_bin2_2,cheb3_bkgPDF_bin2_0,cheb3_bkgPDF_bin2_1,sumexp2_bkgPDF_bin2_x0,sumexp2_bkgPDF_bin2_x1,sumexp2_bkgPDF_bin2_c0,sumexp2_bkgPDF_bin2_c1,pdfindex_bin2,shapeBkg_bkg_bin2__norm,shapeSig_sgn_bin2__norm
    python3 plot_bias_pull.py -i higgsCombine.bias_truth_index_${GEN_PDF_IDX}_fit_envelop_${NAME}.MultiDimFit.mH91.123456.root -o ${NAME}_truth_${GEN_PDF_IDX}_fit_envelope --n-toys $TOYS_FOR_BIAS --r-truth $R_FOR_BIAS --gen-fnc PDF_${GEN_PDF_IDX} --fit-fnc all
   
    mv ${NAME}_truth_${GEN_PDF_IDX}_fit_envelope_pull.png $ORGDIR
    echo "${NAME}_truth_${GEN_PDF_IDX}_fit_envelope_pull.png has been created"
    mv ${NAME}_truth_${GEN_PDF_IDX}_fit_envelope_pull_withfit.png $ORGDIR
    echo "${NAME}_truth_${GEN_PDF_IDX}_fit_envelope_pull_withfit.png has been created"
    

fi


if [[ "${PROCESS}" == "generate" ]]
then

    GEN_PDF_IDX=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;
    
    INDEX_TOY=""
    if [[ "${BRANCH}" == "" ]]
    then
      INDEX_TOY="--setParameters pdfindex_bin1=$GEN_PDF_IDX,pdfindex_bin2=$GEN_PDF_IDX,pdfindex_bin3=$GEN_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else 
       INDEX_TOY="--setParameters $BRANCH=$GEN_PDF_IDX --freezeParameters $BRANCH"
    fi 
 
    combine -M GenerateOnly param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .generate_truth_index_${GEN_PDF_IDX}_${NAME}  --expectSignal $R_FOR_BIAS --saveToys ${INDEX_TOY}
    
    mv higgsCombine.generate_truth_index_${GEN_PDF_IDX}_${NAME}.GenerateOnly.mH91.123456.root $ORGDIR/higgsComb_toys_${GEN_PDF_IDX}_${NAME}.root
    echo "higgsComb_toys_${GEN_PDF_IDX}_${NAME}.root has been created"

fi

if [[ "${PROCESS}" == "bias-envelope-infile" ]]
then
    INPUT=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;
    
    INDEX_FIT=""
    if [[ "${BRANCH}" == "" ]]
    then
       INDEX_FIT="--saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else 
       INDEX_FIT="--saveSpecifiedIndex $BRANCH --cminRunAllDiscreteCombinations"
    fi 

#    echo "Run bias for the envelope (gen pdf from envelope), pdf index for toys = $GEN_PDF_IDX, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"


    cp ${ORGDIR}/${INPUT} .
    combine -M MultiDimFit param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .truth_file_fit_envelop_${NAME} --toysFile ${INPUT} --algo singles  --rMin -100 --rMax 100 --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace ${INDEX_FIT} 
    python3 plot_bias_pull.py -i higgsCombine.truth_file_fit_envelop_${NAME}.MultiDimFit.mH91.123456.root -o ${NAME}_truth_file_fit_envelop --n-toys $TOYS_FOR_BIAS --r-truth $R_FOR_BIAS --gen-fnc External --fit-fnc all
   
    mv ${NAME}_truth_file_fit_envelop_pull.png $ORGDIR
    echo "${NAME}_truth_file_fit_envelop_pull.png has been created"
    mv ${NAME}_truth_file_fit_envelop_pull_withfit.png $ORGDIR
    echo "${NAME}_truth_file_fit_envelop_pull_withfit.png has been created"
    
fi



if [[ "${PROCESS}" == "bias-pdf-infile" ]]
then
    
    INPUT=$1
    shift;
    FIT_PDF_IDX=$1
    shift;
    R_FOR_BIAS=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;

    INDEX_FIT=""

    if [[ "${BRANCH}" == "" ]]
    then
      INDEX_FIT="--setParameters pdfindex_bin1=$FIT_PDF_IDX,pdfindex_bin2=$FIT_PDF_IDX,pdfindex_bin3=$FIT_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3 --saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else 
       INDEX_FIT="--setParameters $BRANCH=$FIT_PDF_IDX --freezeParameters $BRANCH --saveSpecifiedIndex $BRANCH"
    fi 

    #echo "Run bias of idividual padfs on the envelope, pdf index for toys = $GEN_PDF_IDX, pdf index for fit = $FIT_PDF_IDX, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"
    echo $INDEX_FIT

    cp ${ORGDIR}/${INPUT} . 
    combine -M MultiDimFit param_datacard_${NAME}.root --algo singles -m 91 -t $TOYS_FOR_BIAS -n .bias_truth_file_fit_index_${FIT_PDF_IDX}_${NAME} --toysFile ${INPUT} --rMin -100 --rMax 100 --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace --saveToys ${INDEX_FIT}

#--trackParameters cheb3_bkgPDF_bin1_0,cheb3_bkgPDF_bin1_1,cheb3_bkgPDF_bin1_2,cheb3_bkgPDF_bin1_0,cheb3_bkgPDF_bin1_1,sumexp2_bkgPDF_bin1_x0,sumexp2_bkgPDF_bin1_x1,sumexp2_bkgPDF_bin1_c0,sumexp2_bkgPDF_bin1_c1,pdfindex_bin1,shapeBkg_bkg_bin1__norm,shapeSig_sgn_bin1__norm

    python3 plot_bias_pull.py -i higgsCombine.bias_truth_file_fit_index_${FIT_PDF_IDX}_${NAME}.MultiDimFit.mH91.123456.root -o ${NAME}_truth_file_fit_${FIT_PDF_IDX} --n-toys $TOYS_FOR_BIAS --r-truth $R_FOR_BIAS --gen-fnc external --fit-fnc PDF_${FIT_PDF_IDX}
   
    mv ${NAME}_truth_file_fit_${FIT_PDF_IDX}_pull.png $ORGDIR
    echo "${NAME}_truth_file_fit_${FIT_PDF_IDX}_pull.png has been created"
    mv ${NAME}_truth_file_fit_${FIT_PDF_IDX}_pull_withfit.png $ORGDIR
    echo "${NAME}_truth_file_fit_${FIT_PDF_IDX}_pull_withfit.png has been created"
    

fi


if [[ "${PROCESS}" == "bias-scan-envelope" ]]
then
    INPUT=$1
    shift;
    RMAX=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;
    
    INDEX_FIT=""
    if [[ "${BRANCH}" == "" ]]
    then
       INDEX_FIT="--saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else 
       INDEX_FIT="--saveSpecifiedIndex $BRANCH"
    fi 

#    echo "Run bias for the envelope (gen pdf from envelope), pdf index for toys = $GEN_PDF_IDX, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"


    cp ${ORGDIR}/${INPUT} .
    combine -M MultiDimFit param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .bias_grid_truth_file_fit_envelop_${NAME} --toysFile ${INPUT} --algo grid --points 20  --rMin -${RMAX} --rMax ${RMAX} --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace ${INDEX_FIT}  --cminRunAllDiscreteCombinations
    
fi

if [[ "${PROCESS}" == "bias-scan-pdf" ]]
then
    INPUT=$1
    shift;
    FIT_PDF_IDX=$1
    shift;
    RMAX=$1
    shift;
    TOYS_FOR_BIAS=$1
    shift;
    BRANCH=$1
    shift;


    if [[ "${BRANCH}" == "" ]]
    then
      INDEX_FIT="--setParameters pdfindex_bin1=$FIT_PDF_IDX,pdfindex_bin2=$FIT_PDF_IDX,pdfindex_bin3=$FIT_PDF_IDX --freezeParameters pdfindex_bin1,pdfindex_bin2,pdfindex_bin3 --saveSpecifiedIndex pdfindex_bin1,pdfindex_bin2,pdfindex_bin3"
    else
      INDEX_FIT="--setParameters $BRANCH=$FIT_PDF_IDX --freezeParameters $BRANCH --saveSpecifiedIndex $BRANCH"
    fi
    
#    echo "Run bias for the envelope (gen pdf from envelope), pdf index for toys = $GEN_PDF_IDX, r = $R_FOR_BIAS and Ntoys = $TOYS_FOR_BIAS"


    cp ${ORGDIR}/${INPUT} .
    combine -M MultiDimFit param_datacard_${NAME}.root -m 91 -t $TOYS_FOR_BIAS -n .bias_grid_truth_file_fit_index_${FIT_PDF_IDX}_${NAME} --toysFile ${INPUT} --algo grid --points 20  --rMin -${RMAX} --rMax ${RMAX} --saveNLL --cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --saveWorkspace ${INDEX_FIT}
    
fi



rm param_datacard_${NAME}.root

if (( $CLEANCOMB == 1 )) 
then
  echo "cleaning combine space"
  rm workspace_v*.root
  rm param_datacard_*.root
  rm higgsCombine.*.root
fi
cd $ORGDIR

