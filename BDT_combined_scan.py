import ROOT as rt
import os




name_txt = "run2018_bdt_cut_count.txt"
root_name="2018"


os.system("rm -I "+name_txt)
card="imax    1 number of bins\n"
card+="jmax    1 number of processes minus 1\n"
card+="kmax    * number of nuisance parameters\n"
card+="-"*80+"\n"
card+="-"*80+"\n"
card+="bin          signal_region\n"
card+="observation  599\n"
card+="-"*80+"\n"
card+="bin          signal_region  signal_region\n"
card+="process      bkg            signal\n"
card+="process      1              0\n"

current_bf=init_branch_ratio
current_yld=yield_at_upper_init

for i in range(20):
  current_card=card
  current_card+="rate     599            {0}".format(str(current_yld))
  current_bf*=0.9
  current_yld*=0.9
  with open("temp_datacrd.txt",'w') as dc:
     dc.write(current_card)
  os.system("combine -M AsymptoticLimits temp_datacrd.txt -m {0}  -n {1} >> {2}".format(current_bf,root_name,name_txt))
    






