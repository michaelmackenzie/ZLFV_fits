import ROOT as rt
import os

fit_type_1bin=True
init_branch_ratio=7.5e-7
yield_at_upper_init=171
BKG_yield="635"
name_txt = "run2018_shapes.txt"
root_name="2018"

#shape fit only
root_shapes="fit_toys_lepm_background_test.root"
signal_name="$PROCESS_7.5e-7times0.9in$MASS"
bkg_name="$PROCESS"



os.system("rm -I "+name_txt)
card="imax    1 number of bins\n"
card+="jmax    1 number of processes minus 1\n"
card+="kmax    * number of nuisance parameters\n"
card+="-"*80+"\n"
if not fit_type_1bin:
   card+="shapes    hsignal * "+root_shapes+"  "+signal_name+"\n"
   card+="shapes    * * "+root_shapes+"  "+bkg_name+"\n"
card+="-"*80+"\n"
card+="bin          signal_region\n"
if fit_type_1bin:
   card+="observation  "+BKG_yield+"\n"
else:
   card+="observation  -1\n"
card+="-"*80+"\n"
card+="bin          signal_region  signal_region\n"
card+="process      hbackground    hsignal\n"
card+="process      1              0\n"


if fit_type_1bin:
  current_bf=init_branch_ratio
  current_yld=yield_at_upper_init

  for i in range(20):
    current_card=card
    current_card+="rate    "+BKG_yield+"     {0}".format(str(current_yld))
    current_bf*=0.9
    current_yld*=0.9
    with open("temp_datacard.txt",'w') as dc:
       dc.write(current_card)
    os.system("combine -M AsymptoticLimits temp_datacard.txt -m {0}  -n {1} >> {2}".format(current_bf,root_name,name_txt)) 
else:
  card+="rate    -1     -1"
  with open("temp_datacard.txt",'w') as dc:
    dc.write(card)
  for i in range(10):
    os.system("combine -M AsymptoticLimits temp_datacard.txt -m {0}  -n {1} >> {2}".format(str(i),root_name,name_txt))






