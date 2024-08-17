# LFV Z resonance fits

## Z prime scan

The Z prime scan searches for a narrow resonance in the e-mu data using the Z->e+mu analysis framework and BDT.

### Create signal and background PDFs for a range of mass points for a single BDT category
```
MINMASS=110
MAXMASS=500
MINBDT=0.70
MAXBDT=1.01
NAME=bdt_0d7_1d0_v01
time python ScanMuE_fit_wrapper_v1.py -o ${NAME} --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min ${MINBDT} --xgb-max ${MAXBDT}
ls -l datacards/${NAME}/combine_combine_zprime_${NAME}_mp*.txt | head -n 2
ls -l WorkspaceScanSGN/workspace_scansgn_v2_${NAME}_mp*.root | head -n 2
ls -l WorkspaceScanBKG/workspace_scanbkg_v2_${NAME}_mp*.root | head -n 2
#figures are printed to: figures/${NAME}/ (signal) and figures/${NAME}_mp*/ (background/data)
```

### Create standard BDT categories
```
MINMASS=110
MAXMASS=500
NAME=v01
time ./make_scan_cards.sh --min-mass ${MINMASS} --max-mass ${MAXMASS} --tag ${NAME}
ls -l datacards/bdt_${NAME}/combine_combine_zprime_${NAME}_mp*.txt | head -n 2
```

### Scan the mass points, evaluating signal rates and upper limits
```
NAME=v01
time python perform_scan.py -o bdt_${NAME} [--asimov] [--unblind]
ls -l figures/scan_bdt_${NAME}[_asimov]/*.png
```
