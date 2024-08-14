import os
import argparse
import ROOT as rt
from array import array
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')

def plot_graph(mpoint_array,signal_array,name,xaxis="",yaxis=""):
   signal_vs_mass = rt.TGraph(len(mpoint_array), mpoint_array, signal_array)

   fnc = rt.TF1("fnc", "pol2", 200, 1000)
   par_yld = array('d',[0, 0, 0])
   cnv = rt.TCanvas("cnv_"+name,"",800,600)
   signal_vs_mass.Draw("A*")
   signal_vs_mass.SetTitle("#bf{CMS} #it{Preliminary};"+xaxis+";"+yaxis)
   signal_vs_mass.Fit(fnc,"R")
   cnv.SaveAs(name+".png")
   fnc.GetParameters(par_yld)
   return par_yld


parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="test", type=str,help="output root name")
parser.add_argument("--data-file", dest="data_file",default="", type=str,help="data file")
parser.add_argument("--xgb-min", dest="xgb_min",default="0.70", type=str,help="data file")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.01", type=str,help="data file")
parser.add_argument("--scan-min", dest="scan_min",default=110., type=float,help="mass signal")
parser.add_argument("--scan-max", dest="scan_max",default=500., type=float,help="mass signal")
parser.add_argument("--scan-step", dest="scan_step",default=0.5, type=float,help="mass signal")
parser.add_argument("--component", dest="component",type=str,default="all",help="Options: all sgn bkg ")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="data file")
parser.add_argument("--skip-shape-dc", dest="skip_shape_dc",default=False, action='store_true',help="shape experiment")
parser.add_argument("--year", dest="year",default="Run2", type=str,help="year")
parser.add_argument("--fit-version", dest="ver",default="2", type=str,help="fit version")
parser.add_argument("--outvar", dest="outvar",default="mass_ll",type=str,help="data file")
parser.add_argument("--param-name", dest="param_name",default="bin",type=str,help="data file")
parser.add_argument("--skip-fit", dest="skip_fit",default=False, action='store_true',help="fit skip")
parser.add_argument("--skip-sgn-syst", dest="skip_sgn_syst",default=False, action='store_true',help="shape experiment")
parser.add_argument("--skip-bkg-altfits", dest="skip_bkg_altfits",default=False, action='store_true',help="shape experiment")

args, unknown = parser.parse_known_args()


##### configuaration ########

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exitting"
   exit()

MaxMasspoints=-1 #X for debug; -1 for run

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"

### MC signal mass points
sgn_masspoints=["200","400","600","800","1000"]
sgn_masspoint_files={
               "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
               "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
               "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}
ndens={"200":44844.,"400":55294.,"600":60192.,"800":64756.,"1000":67263.}

### Data
data_files={"2016":"Meas_full_bdt_v7_emu_scan_data_Run16.root","2017":"Meas_full_bdt_v7_emu_scan_data_Run17.root","2018":"Meas_full_bdt_v7_emu_scan_data_Run18.root","Run2":"Meas_full_bdt_v7_emu_scan_data_Run1*.root"}
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}

### sf
#sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)"
sf="1"


### arange flags
unblind="false"
if args.unblind:
   unblind="true"
shape_dc="true"
if args.skip_shape_dc:
   shape_dc="false"
do_sgn_syst="true"
if args.skip_sgn_syst:
   do_sgn_syst="false"
do_bkg_altfits="true"
if args.skip_bkg_altfits:
   do_bkg_altfits="false"
if args.data_file=="":
   args.data_file=path+data_files[args.year]

rt.gROOT.SetBatch(True)



##### mass points analysis ######
### efficiency and yield vs mass
nsgn_array = array('d',[])
toteff_array=  array('d',[])
bdteff_array=  array('d',[])
mpoint_array = array('d',[])
for mpoint in sgn_masspoints:
  cc=rt.TChain("mytreefit")
  cc.Add(path+"/"+sgn_masspoint_files[mpoint])
  htemp = rt.TH1F("htemp_"+mpoint,"",1,0,2)
  cc.Draw("1>>htemp_"+mpoint,sf+"*(Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<=xgb && xgb<"+args.xgb_max+")")
  nnum = htemp.Integral()
  nsgn_array.append(lumis[args.year]*1*nnum/ndens[mpoint])
  mpoint_array.append(float(mpoint))
  toteff_array.append(nnum/(1.0*ndens[mpoint]))
  bdteff_array.append(nnum/(1.0*cc.GetEntries()))
  print "total pass ("+mpoint+") =",nnum,"eff",nnum/(1.0*ndens[mpoint]),"yield",nsgn_array[-1]

par_yield = plot_graph(mpoint_array,nsgn_array,"yld_vs_mass","m(e,#mu)","Yield")
par_bdt_eff = plot_graph(mpoint_array,bdteff_array,"bdteff_vs_mass","m(e,#mu)","BDT eff.")
par_total_eff = plot_graph(mpoint_array,toteff_array,"totaleff_vs_mass","m(e,#mu)","Total eff.")

### get widths
width_array = array('d',[])
for mpoint in sgn_masspoints:
  cc=rt.TChain("mytreefit")
  cc.Add(path+"/"+sgn_masspoint_files[mpoint])
  htemp = rt.TH1F("htemp_"+mpoint,"",50,float(mpoint)*0.8,float(mpoint)*1.2)
  cc.Draw("mass_ll>>htemp_"+mpoint,sf+"*(Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<=xgb && xgb<"+args.xgb_max+")")

  gauss = rt.TF1("gauss", "gaus", float(mpoint)*0.95,float(mpoint)*1.05)
  par_gs = array('d',[0, 0, 0])
  cnv = rt.TCanvas("cnv_"+mpoint,"",800,600)
  htemp.Draw("HIST")
  htemp.Fit(gauss,"LR")
  htemp.SetTitle("#bf{CMS} #it{Preliminary};m(e,#mu);")
  gauss.Draw("sames")
  cnv.SaveAs("mass_fit_"+mpoint+".png")
  gauss.GetParameters(par_gs)
  width_array.append(par_gs[2])
  
par_width = plot_graph(mpoint_array,width_array,"width_vs_mass","m(e,#mu)","Width")

###### main code ######
if MaxMasspoints>0:
   print "WILL PRODUCE",MaxMasspoints,"points"
NextPoint=True 

# set intitial scan mass
sr_center=args.scan_min

cnt=-1
while (NextPoint):
  cnt+=1

  # calculate blind range and expected widths/yields for a mass
  sr_width = round(par_width[0]+par_width[1]*sr_center+par_width[2]*sr_center*sr_center,2)
  sr_min = round(sr_center - 1*sr_width,2)
  sr_max = round(sr_center + 1*sr_width,2)
  sr_yld = round(par_yield[0] + par_yield[1]*sr_center + par_yield[2]*sr_center*sr_center,2)
  min_mass = round(sr_center - 5*sr_width,2)
  max_mass = round(sr_center + 5*sr_width,2)
  print "SR central",sr_center,"width",sr_width,"min",sr_min,"max",sr_max,"yield",sr_yld

  # create pdfs for mass point
  if not args.skip_fit:
     if args.component == "sgn" or args.component == "all":
        os.system('root -l -b -q ScanMuE_fit_sgn_v'+args.ver+'.C\'("'+args.name+"_mp"+str(cnt)+'",'+str(min_mass)+','+str(max_mass)+','+str(sr_center)+','+str(sr_width)+','+str(sr_yld)+','+shape_dc+',"'+args.outvar+'",'+do_sgn_syst+',"'+args.param_name+'")\'')
     if args.component == "bkg" or args.component == "all":
        os.system('root -l -b -q ScanMuE_fit_bkg_v'+args.ver+'.C\'("'+args.name+"_mp"+str(cnt)+'","'+args.data_file+'","'+args.xgb_min+'","'+args.xgb_max+'",'+str(min_mass)+','+str(max_mass)+','+str(sr_min)+','+str(sr_max)+','+unblind+','+shape_dc+',"'+args.outvar+'","'+args.param_name+'")\'')

  # next iteration mass and exit conditions
  sr_center = round(sr_center +args.scan_step*sr_width,2)
  if cnt>= MaxMasspoints and MaxMasspoints>0:
     print "Requested only",MaxMasspoints,"run"
     NextPoint=False
  if sr_center >= args.scan_max:
     NextPoint=False
  
  
print "scannned",cnt,"points"


