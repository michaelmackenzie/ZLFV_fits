# Evaluate systematic uncertainties for the Z prime signal
import os
import argparse
import ROOT as rt
from array import array

def eval_unc(tree, name, cuts, mass):
    c = rt.TCanvas()
    h_nom = rt.TH1F("h_nom","",1000,mass-100,mass+100)
    cc.Draw("mass_ll>>h_nom",cuts+"*"+name+"_wt")
    h_up = rt.TH1F("h_up","",1000,mass-100,mass+100)
    cc.Draw("mass_ll>>h_up",cuts+"*"+name+"_up")
    h_down = rt.TH1F("h_down","",1000,mass-100,mass+100)
    cc.Draw("mass_ll>>h_down",cuts+"*"+name+"_down")

    nom  = h_nom.Integral()
    up   = h_up.Integral()
    down = h_down.Integral()

    if nom <= 0.:
        cc.Draw("mass_ll>>h_nom",cuts)
        nom  = h_nom.Integral()
        
    # See if it's defined by 'name_sys' instead
    if up <= 0. or down <= 0.:
        cc.Draw("mass_ll>>h_up",cuts+"*"+name+"_sys")
        up = h_up.Integral()
        if up > 0.:
            down = (nom/up) * nom
        

    return [nom, up, down]

rt.gROOT.SetBatch(True)

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--xgb-min", dest="xgb_min",default="0.70", type=str,help="data file")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.01", type=str,help="data file")
parser.add_argument("--mass-point", dest="mass_point",default="", type=str,help="Mass point to process")

args, unknown = parser.parse_known_args()
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()


### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/signal_unc/"
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

### MC signal mass points
sgn_masspoints=["200","400","600","800","1000"]
sgn_masspoint_files={
               "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
               "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
               "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}
ndens={"200":96300.,"400":97600.,"600":97600.,"800":97600.,"1000":97600.}

# Define the uncertainties
uncs = ['SFbtag', 'Electron_dxydz', 'Muon_dxydz', 'Trg', 'Prefire', 'PU', 'Electron_IsoID', 'Muon_IsoID',
        'Electron_ID', 'Muon_ID', 'Electron_RecoID', 'Muon_RecoID',
        'PtZ', 'PDF', 'Scale']

# Selection cuts
cuts = "(Flag_met && Flag_muon && Flag_electron) && xgb >= %s && xgb < %s" % (args.xgb_min, args.xgb_max)

# Loop through the signal files
effects = []
for mpoint in sgn_masspoints:
    if args.mass_point != "" and args.mass_point != mpoint: continue
    print "Processing mass point", mpoint
    # Create a TChain for the signal
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_masspoint_files[mpoint])
    # Loop through each uncertainty
    for unc in uncs:
        [nom, up, down] = eval_unc(cc, unc, cuts, int(mpoint))
        effects.append([up/nom, down/nom])

        
    for index in range(len(uncs)):
        print "%-15s: %.3f/%.3f" % (uncs[index], effects[index][0], effects[index][1])

