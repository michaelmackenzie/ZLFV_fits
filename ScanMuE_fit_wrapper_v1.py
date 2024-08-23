# Perform a scan of the emu data, fitting for a potential signal resonance
import os
import argparse
import ROOT as rt
from array import array
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')

#----------------------------------------------------------------------------------------
# Write a combine data card for a single bin
def print_datacard(name, sig_file, bkg_file, param_name, mass):
   # Re-create the card
   f = open(name, "w")

   # Standard preamble
   f.write("# -*- mode: tcl -*-\n")
   f.write("#Auto generated Z prime search COMBINE datacard\n")
   f.write("#Using Z prime mass = %.2f\n\n" % (mass))
   f.write("imax 1 number of channels\njmax * number of backgrounds\nkmax * number of nuisance parameters\n\n")

   # Define the background and signal PDFs
   f.write("-----------------------------------------------------------------------------------------------------------\n")
   f.write("shapes signal * %s workspace_signal:signal_pdf_%s\n" % (sig_file, param_name))
   f.write("shapes background * %s ws_bkg:multipdf_%s\n" % (bkg_file, param_name))
   f.write("shapes data_obs * %s ws_bkg:data_obs\n" % (bkg_file))
   f.write("-----------------------------------------------------------------------------------------------------------\n\n")

   # Define the channel
   f.write("bin         bin\n")
   f.write("observation  -1\n\n")
   f.write("bin         bin       bin\n")
   f.write("process    signal   background\n\n")
   f.write("process       0         1\n\n")
   f.write("rate          1         1\n") #Rate is taken from the _norm variables

   # Define the uncertainties
   f.write("-----------------------------------------------------------------------------------------------------------\n")
   #FIXME: Implement mass-dependent uncertainties
   f.write("ElectronID lnN 1.02     -\n")
   f.write("MuonID     lnN 1.02     -\n")
   f.write("Lumi       lnN 1.02     -\n")
   f.write("BTag       lnN 1.005    -\n")
   f.write("Theory     lnN 1.01     -\n")
   f.write("BDT        lnN 1.02     -\n")
   f.write("-----------------------------------------------------------------------------------------------------------\n\n")

   # Scale uncertainties
   f.write("-----------------------------------------------------------------------------------------------------------\n")
   f.write("elec_ES_shift_%s param 0 1 [-7, 7]\n" % (param_name))
   f.write("muon_ES_shift_%s param 0 1 [-7, 7]\n" % (param_name))
   f.write("-----------------------------------------------------------------------------------------------------------\n\n")

   # Define the envelope discrete index to be scanned
   f.write("pdfindex_%s discrete\n" % (param_name))

   # Close the file
   f.close()
   

#----------------------------------------------------------------------------------------
# Make a plot and save the figure, fit a 2nd order polynomial to it
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


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

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
parser.add_argument("--skip-dc", dest="skip_dc",default=False, action='store_true',help="Skip datacard creation")
parser.add_argument("--mass-point", dest="mass_point",default=-1,type=int,help="Single mass point to process")
parser.add_argument("--full-mass", dest="full_mass",default=False,action='store_true',help="Fit the entire mass distribution")
parser.add_argument("--log-files", dest="log_files",default=False,action='store_true',help="Write individual mass point fits to log files")
parser.add_argument("--dry-run", dest="dry_run",default=False,action='store_true',help="Don't execute script calls")

args, unknown = parser.parse_known_args()


##### configuaration ########

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exitting"
   exit()

MaxMasspoints=-1 #X for debug; -1 for run

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/%s/" % (args.name)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))
os.system("[ ! -d log ] && mkdir log")
if not os.path.exists(carddir+"WorkspaceScanBKG"):
   os.symlink("../../WorkspaceScanBKG", carddir+"WorkspaceScanBKG")
if not os.path.exists(carddir+"WorkspaceScanSGN"):
   os.symlink("../../WorkspaceScanSGN", carddir+"WorkspaceScanSGN")
   

### MC signal mass points
sgn_masspoints=["200","400","600","800","1000"]
sgn_masspoint_files={
               "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
               "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
               "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}
ndens={"200":96300.,"400":55294.,"600":60192.,"800":64756.,"1000":67263.}

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


#----------------------------------------------
# Get the signal parameterization vs. mass
#----------------------------------------------

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
  cross_section = 1. #in units femto-barns * BR(Z'->emu)
  nsgn_array.append(lumis[args.year]*1*nnum/ndens[mpoint])
  mpoint_array.append(float(mpoint))
  toteff_array.append(nnum/(1.0*ndens[mpoint]))
  bdteff_array.append(nnum/(1.0*cc.GetEntries()))
  print "total pass ("+mpoint+") =",nnum,"eff",nnum/(1.0*ndens[mpoint]),"yield",nsgn_array[-1]

par_yield = plot_graph(mpoint_array,nsgn_array,figdir+"yld_vs_mass","m(e,#mu)","Yield")
par_bdt_eff = plot_graph(mpoint_array,bdteff_array,figdir+"bdteff_vs_mass","m(e,#mu)","BDT eff.")
par_total_eff = plot_graph(mpoint_array,toteff_array,figdir+"totaleff_vs_mass","m(e,#mu)","Total eff.")

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
  cnv.SaveAs(figdir+"mass_fit_"+mpoint+".png")
  gauss.GetParameters(par_gs)
  width_array.append(par_gs[2])
  
par_width = plot_graph(mpoint_array,width_array,figdir+"width_vs_mass","m(e,#mu)","Width")

#----------------------------------------------
# Perform the signal scan at all mass points
#----------------------------------------------

###### main code ######
if MaxMasspoints>0:
   print "WILL PRODUCE",MaxMasspoints,"points"
NextPoint=True 

# set intitial scan mass
sr_center=args.scan_min

cnt=-1
script_head = 'echo ' if args.dry_run else ''
while (NextPoint):
  cnt+=1

  # calculate blind range and expected widths/yields for a mass
  sr_buffer = 1 #in signal width units
  region_buffer = 10 #in signal width units
  sr_width = round(par_width[0]+par_width[1]*sr_center+par_width[2]*sr_center*sr_center,2)
  sr_min = round(sr_center - sr_buffer*sr_width,2)
  sr_max = round(sr_center + sr_buffer*sr_width,2)
  sr_yld = round(par_yield[0] + par_yield[1]*sr_center + par_yield[2]*sr_center*sr_center,2)
  min_mass = round(sr_center - region_buffer*sr_width,2) if not args.full_mass else  90.
  max_mass = round(sr_center + region_buffer*sr_width,2) if not args.full_mass else 700.
  print "SR central",sr_center,"width",sr_width,"min",sr_min,"max",sr_max,"yield",sr_yld

  if(args.mass_point < 0 or cnt == args.mass_point):
    # create pdfs for mass point
    if not args.skip_fit:
      if args.component == "sgn" or args.component == "all":
        tail = (' >| log/fit_sgn_%s_mp%i.log' % (args.name, cnt)) if args.log_files else ''
        os.system(script_head + 'root -l -b -q ScanMuE_fit_sgn_v'+args.ver+'.C\'("'
                  +args.name+"_mp"+str(cnt)+'",'
                  +str(min_mass)+','+str(max_mass)+','+str(sr_center)+','+str(sr_width)+','
                  +str(sr_yld)+','+shape_dc+',"'+args.outvar+'",'+do_sgn_syst+',"'+args.param_name+'")\'' + tail)
      if args.component == "bkg" or args.component == "all":
         tail = (' >| log/fit_bkg_%s_mp%i.log' % (args.name, cnt)) if args.log_files else ''
         os.system(script_head + 'root -l -b -q ScanMuE_fit_bkg_v'+args.ver+'.C\'("'+args.name+"_mp"+str(cnt)+'","'+args.data_file
                   +'","'+args.xgb_min+'","'+args.xgb_max+'",'+str(min_mass)+','+str(max_mass)+','+str(sr_min)+','+str(sr_max)
                   +','+unblind+','+shape_dc+',"'+args.outvar+'","'+args.param_name+'")\'' + tail)
        

    # Create a corresponding datacard
    cardname = carddir + "combine_zprime_" + args.name + "_mp" + str(cnt) + ".txt"
    sig_file = "WorkspaceScanSGN/workspace_scansgn_v" + args.ver + "_" + args.name + "_mp" + str(cnt) + ".root"
    bkg_file = "WorkspaceScanBKG/workspace_scanbkg_v" + args.ver + "_" + args.name + "_mp" + str(cnt) + ".root"
    if not args.dry_run: print_datacard(cardname, sig_file, bkg_file, args.param_name, sr_center)


  # next iteration mass and exit conditions
  sr_approx_width = sr_center*(4./200.) # approximate as a linear function so it's stable between BDT categories where the width may vary
  sr_center = round(sr_center +args.scan_step*sr_approx_width,2)
  if cnt>= MaxMasspoints and MaxMasspoints>0:
     print "Requested only",MaxMasspoints,"run"
     NextPoint=False
  if sr_center >= args.scan_max:
     NextPoint=False
  
  
print "scannned",cnt,"points"


