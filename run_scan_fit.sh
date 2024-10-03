MINMASS=110
MAXMASS=500
COMPONENT=all 


time python ScanMuE_fit_wrapper_v2.py -o bdt_bin1_v2 --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min 0.3 --xgb-max 0.7 --component ${COMPONENT} --outvar lepm_11 --param-name bin1
time python ScanMuE_fit_wrapper_v2.py -o bdt_bin2_v2 --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min 0.7 --component ${COMPONENT} --outvar lepm_12 --param-name bin2

#time python ScanMuE_fit_wrapper_v1_michl.py -o bdt_bin1_v2 --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min 0.3 --xgb-max 0.7 --component ${COMPONENT} --outvar lepm_11 --param-name bin1 --skip-fit
#time python ScanMuE_fit_wrapper_v1_michl.py -o bdt_bin2_v2 --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min 0.7 --component ${COMPONENT} --outvar lepm_12 --param-name bin2 --skip-fit


ls -l DatacardsScan/${NAME}/combine_combine_zprime_${NAME}_mp*.txt | head -n 2
ls -l WorkspaceScanSGN/workspace_scansgn_v2_${NAME}_mp*.root | head -n 2
ls -l WorkspaceScanBKG/workspace_scanbkg_v2_${NAME}_mp*.root | head -n 2
#figures are printed to: figures/${NAME}/ (signal) and figures/${NAME}_mp*/ (background/data)
