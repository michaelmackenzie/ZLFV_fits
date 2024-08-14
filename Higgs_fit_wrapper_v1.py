import os
import argparse
import ROOT as rt
from array import array
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')


parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="test", type=str,help="output root name")
parser.add_argument("--data-file", dest="data_file",default="", type=str,help="data file")
parser.add_argument("--xgb-min", dest="xgb_min",default="0.75", type=str,help="data file")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.01", type=str,help="data file")
parser.add_argument("--mass", dest="mass",default=125., type=float,help="mass signal")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="data file")
parser.add_argument("--create-shape-dc", dest="shape_dc",default=False, action='store_true',help="shape experiment")
parser.add_argument("--lumi", dest="lumi",default=0.0, type=float,help="lumi to compute signal")
parser.add_argument("--year", dest="year",default="Run2", type=str,help="year")
parser.add_argument("--fit-version", dest="ver",default="1", type=str,help="fit version")
parser.add_argument("--outvar", dest="outvar",default="mass_ll",type=str,help="data file")
parser.add_argument("--param-name", dest="param_name",default="bin",type=str,help="data file")
parser.add_argument("--skip-fit", dest="skip_fit",default=False, action='store_true',help="fit skip")
parser.add_argument("--skip-sgn-syst", dest="skip_sgn_syst",default=False, action='store_true',help="shape experiment")
parser.add_argument("--skip-bkg-altfits", dest="skip_bkg_altfits",default=False, action='store_true',help="shape experiment")

args, unknown = parser.parse_known_args()

if len(unknown)>0: 
   print "not found:",unknown,"exitting"
   exit()

path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"

sgn_masspoint_file="Meas_full_bdt_v7_higgs_emu_mcRun18.root"
nden=40000.

unblind="false"
if args.unblind:
   unblind="true"
shape_dc="false"
if args.shape_dc:
   shape_dc="true"
do_sgn_syst="true"
if args.skip_sgn_syst:
   do_sgn_syst="false"
do_bkg_altfits="true"
if args.skip_bkg_altfits:
   do_bkg_altfits="false"


data_files={"2016":"Meas_full_bdt_v7_zprime_data_emu_Run16.root","2017":"Meas_full_bdt_v7_zprime_data_emu_Run17.root","2018":"Meas_full_bdt_v7_zprime_data_emu_Run18.root","Run2":"Meas_full_bdt_v7_zprime_data_emu_Run1*.root"}
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}

#sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)"
sf="1"

if args.lumi==0:
   args.lumi=lumis[args.year]
if args.data_file=="":
   args.data_file=path+data_files[args.year]

rt.gROOT.SetBatch(True)

nsgn=0.
toteff=0.
bdteff=0.
cc=rt.TChain("mytreefit")
cc.Add(path+"/"+sgn_masspoint_file)
htemp = rt.TH1F("htemp","",1,0,2)
cc.Draw("1>>htemp",sf+"*(Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<=xgb && xgb<"+args.xgb_max+")")
nnum = htemp.Integral()
print "total pass =",nnum
nsgn = lumis[args.year]*(48.61+3.766+0.5071+1.358+0.880)*1000*4.7e-5*nnum/nden
toteff= nnum/(1.0*nden)
bdteff= nnum/(1.0*cc.GetEntries())

print "Higgs result: bdt-min",args.xgb_min,"bdt-max",args.xgb_max,"yield",nsgn,"BDT eff",bdteff,"total eff",toteff


if args.skip_fit:
   exit()


os.system('root -l -b -q Zprime_fit_v'+args.ver+'.C\'("'+args.name+'","'+args.data_file+'",'+str(args.mass)+',"'+args.xgb_min+'","'+args.xgb_max+'",'+unblind+','+shape_dc+','+str(nsgn)+',"'+args.outvar+'",'+do_sgn_syst+','+do_bkg_altfits+',"'+args.param_name+'")\'')
