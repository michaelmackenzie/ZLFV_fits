import os
import argparse



parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="test", type=str,help="output root name")
parser.add_argument("--lepton", dest="lepton",default="lep", type=str,help="ele mu")
parser.add_argument("--data-files", dest="data_file",default=None, type=str,help="data file to run on multiple files needs \* instead of *")
parser.add_argument("--save-histo", dest="save_histo",default="true", type=str,help="IF run on ntuple -> stores the histo in fdata_from_fit.root; options: true, false")
parser.add_argument("--save-shape", dest="save_shape",default=False, action='store_true',help="shape experiment")
parser.add_argument("--cpp-code", dest="cppcode",default="ZLL_fit_v7", type=str,help="c++ code to wrap - Extension .C added automaticaly")
parser.add_argument("--use-histo", dest="histo_path",default="None", type=str,help="Directly using histo and skips ntuple reading; path of histo file")
parser.add_argument("--use-histo-name", dest="histo_name",default="hdata", type=str,help="Directly using histo and skips ntuple reading; name of histo")
args, unknown = parser.parse_known_args()



if args.lepton!="ele" and args.lepton!="mu":
   print ("options: mu (muons) or ele (electrons)")
   exit()

if args.data_file==None and args.lepton=="mu":
   args.data_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v5_tuples/BDT_outputs_v5/Zmumu_SingleMuData18/Meas_bdt_v5_Zmumu_Data_Run18*.root"

if args.data_file==None and args.lepton=="ele":
   args.data_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v5_tuples/BDT_outputs_v5/Zee_SingleEleData18/Meas_bdt_v5_Zee_Data_Run18*.root"



save_shape = "false"
if args.save_shape:
   save_shape="true"

print ("cpp code",args.cppcode)
os.system('root -l -b -q '+args.cppcode+'.C\'("'+args.name+'","'+args.lepton+'","'+args.data_file+'",'+save_shape+','+args.save_histo+',"'+args.histo_path+'","'+args.histo_name+'")\'')


