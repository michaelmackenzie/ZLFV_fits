
####### mark 2
### signal
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component sgn -o bin1 --param-name bin1 --outvar lepm_11 --xgb-min 0.3 --xgb-max 0.7 --create-shape-dc
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component sgn -o bin2 --param-name bin2 --outvar lepm_12 --xgb-min 0.7 --xgb-max 0.9 --create-shape-dc
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component sgn -o bin3 --param-name bin3 --outvar lepm_13 --xgb-min 0.9 --create-shape-dc



### Zmm
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component zmm --zmm-file template_zemu_embed_v3_bin1.root -o bin1 --param-name bin1 --outvar lepm_11 --xgb-min 0.3 --xgb-max 0.7 --create-shape-dc --run-histo 
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component zmm --zmm-file template_zemu_embed_v3_bin2.root -o bin2 --param-name bin2 --outvar lepm_12 --xgb-min 0.7 --xgb-max 0.9 --create-shape-dc --run-histo 
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component zmm --zmm-file template_zemu_embed_v3_bin3.root -o bin3 --param-name bin3 --outvar lepm_13 --xgb-min 0.9 --create-shape-dc --run-histo 


### continous bkg
python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component bkg --bkg-file template_zemu_embed_v3_bin1.root -o bin1_new --param-name bin1 --outvar lepm_11 --xgb-min 0.3 --xgb-max 0.7 --create-shape-dc --run-histo --histo-toy --add-pol-order 1 --add-exp-order 1 --add-plaw-order 1
python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component bkg --bkg-file template_zemu_embed_v3_bin2.root -o bin2_new --param-name bin2 --outvar lepm_12 --xgb-min 0.7 --xgb-max 0.9 --create-shape-dc --run-histo --histo-toy  --add-pol-order 1 --add-exp-order 1 --add-plaw-order 1
python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component bkg --bkg-file template_zemu_embed_v3_bin3.root -o bin3_new --param-name bin3 --outvar lepm_13 --xgb-min 0.9 --create-shape-dc --run-histo --histo-toy  --add-pol-order 1 --add-exp-order 1 --add-plaw-order 1



### data and template
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component data --bkg-file template_zemu_embed_v3_bin1.root -o bin1 --param-name bin1 --outvar lepm_11 --xgb-min 0.3 --xgb-max 0.7 --run-histo
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component data --bkg-file template_zemu_embed_v3_bin2.root -o bin2 --param-name bin2 --outvar lepm_12 --xgb-min 0.7 --xgb-max 0.9 --run-histo
#python ZMuE_fit_mk2_wrapper_v1.py --fit-version 1 --component data --bkg-file template_zemu_embed_v3_bin3.root -o bin3 --param-name bin3 --outvar lepm_13 --xgb-min 0.9 --run-histo

