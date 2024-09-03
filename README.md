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
time python perform_scan.py -o bdt_${NAME} [--asimov] [--unblind] [--smooth-expected]
ls -l figures/scan_bdt_${NAME}[_asimov]/*.png
```

### Generate a toy dataset and run a scan over the toy dataset
This assumes the nominal scan is already processed on data with corresponding COMBINE cards available.

```
# Fit the data in the entire mass range for toy generation (only needed once, all toys can be generated from this initial fit)
python ScanMuE_fit_wrapper_v1.py -o bdt_0d3_0d7_LEE --full-mass --scan-min 300 --scan-max 300.1 --scan-step 1 --xgb-min 0.30 --xgb-max 0.70 --param-name bin1 --component bkg
python ScanMuE_fit_wrapper_v1.py -o bdt_0d7_1d0_LEE --full-mass --scan-min 300 --scan-max 300.1 --scan-step 1 --xgb-min 0.70 --xgb-max 1.01 --param-name bin2 --component bkg

# Generate a single toy dataset for each BDT region
python create_toy.py --fit-file WorkspaceScanBKG/workspace_scanbkg_v2_bdt_0d3_0d7_LEE_mp0.root -o toy_0d3_0d7 --toy 2 --param bin1 --seed 90
python create_toy.py --fit-file WorkspaceScanBKG/workspace_scanbkg_v2_bdt_0d7_1d0_LEE_mp0.root -o toy_0d7_1d0 --toy 2 --param bin2 --seed 90

# Create toy scan COMBINE cards from an existing scan dataset
CARDDIR="datacards/bdt_v03_step_1d0/" #Existing datacards
TOYDIR="datacards/bdt_v03_step_1d0_toy_2/" #Output toy datacard directory
TOYFILE1="WorkspaceScanTOY/toy_file_toy_0d3_0d7_2.root" #Low score toy file
TOYFILE2="WorkspaceScanTOY/toy_file_toy_0d7_1d0_2.root"
./clone_cards_for_toy.sh ${CARDDIR} ${TOYDIR} ${TOYFILE1} ${TOYFILE2}

# Run the scan over the toy datacards
time python perform_scan.py -o bdt_v03_step_1d0_toy_2 --unblind --smooth-expected
```

### Validation studies

Validation studies are performed using the [perform_scan_validation.py](perform_scan_validation.py) tool.

```
# Perform a bias test for each mass point in the scan:
NTOYS=1000
NAME="bdt_v06_step_1d0"
time python perform_scan_validation.py -o ${NAME} -t ${NTOYS} --test bias
ls -l figures/val_${NAME}_bias/pulls.png
ls -l figures/val_${NAME}_bias/${NAME}_mp*_bias.png | head -n 5
```

Underlying tools to perform a single validation check:
- [bemu_bias.sh](tests/bemu_bias.sh): Perform a self-closure bias test on a single data card.
- [bemu_gen_fit_test.sh](tests/bemu_gen_fit_test.sh): Perform a closure test using (potentially) different generation and fit data cards.
- [do_goodness_of_fit.sh](tests/do_goodness_of_fit.sh): Perform a goodness-of-fit test.
- [impacts.sh](tests/impacts.sh): Evaluate impacts (including condor configurations if requested).

### Additional useful tools/studies

- [plot_signal_model.py](tools/plot_signal_model.py): Plot the signal model interpolation, including the line shape and the overall efficiency.
- [plot_zprime_bdt.py](tools/plot_zprime_bdt.py): Plot the Z prime BDT score distribution and CDF as a function of mass.
- [eval_zprime_unc.py](tools/eval_zprime_unc.py): Evaluate standard uncertainties on the Z prime signal efficiency.
- [get_ngen.py](tools/get_ngen.py): Retrieve the signal ntuple normalization info.
- [compare_scans.py](tools/compare_scans.py): Plot the ratio of the expected limit from two (identical mass point) scans as a function of mass to compare their sensitivity.
- [optimize_binning.py](tools/optimize_binning.py): Rough optimization of a single BDT bin using cut-and-count limits without systematic uncertainties.
