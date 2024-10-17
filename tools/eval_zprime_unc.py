# Evaluate systematic uncertainties for the Z prime signal
import os
import argparse
import ROOT as rt
from array import array

def eval_unc(tree, name, cuts, mass):
    c = rt.TCanvas()
    h_nom = rt.TH1F("h_nom","",1000,mass-100,mass+100)
    tree.Draw("mass_ll>>h_nom",cuts+"*"+name+"_wt")
    h_up = rt.TH1F("h_up","",1000,mass-100,mass+100)
    tree.Draw("mass_ll>>h_up",cuts+"*"+name+"_up")
    h_down = rt.TH1F("h_down","",1000,mass-100,mass+100)
    tree.Draw("mass_ll>>h_down",cuts+"*"+name+"_down")

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

def eval_scale_unc(tree, tree_up, tree_down, name, cuts, mass, figdir):
    c = rt.TCanvas()
    h_nom = rt.TH1F("h_nom","",1000,mass-100,mass+100)
    tree.Draw("mass_ll>>h_nom",cuts)
    h_up = rt.TH1F("h_up","",1000,mass-100,mass+100)
    tree_up.Draw("mass_ll>>h_up",cuts)
    h_down = rt.TH1F("h_down","",1000,mass-100,mass+100)
    tree_down.Draw("mass_ll>>h_down",cuts)

    nom  = h_nom.Integral()
    up   = h_up.Integral()
    down = h_down.Integral()

    if nom <= 0. or up <= 0. or down <= 0.:
        print 'Failed to process %s uncertainty: nom = %f, up = %f, down = %f' % (name, nom, up, down)
        return [1., 1., 1.]

    nom_mean  = h_nom .GetMean()
    up_mean   = h_up  .GetMean()
    down_mean = h_down.GetMean()
    print name + " mean effect: up-shift = %.5f, down-shift = %.5f" % ((up_mean - nom_mean)/mass, (down_mean - nom_mean)/mass)

    c = rt.TCanvas()
    h_nom.Draw('hist')
    h_up.SetLineColor(rt.kRed)
    h_up.Draw('E1 same')
    h_down.SetLineColor(rt.kGreen)
    h_down.Draw('E1 same')
    h_nom.GetXaxis().SetRangeUser(mass - 0.10*mass, mass + 0.10*mass)
    c.SaveAs(figdir + 'sys_' + name + '_mass_' + str(mass) + '.png')
    return [nom, up, down]

rt.gROOT.SetBatch(True)

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--xgb-min", dest="xgb_min",default="0.70", type=str,help="data file")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.01", type=str,help="data file")
parser.add_argument("--mass-point", dest="mass_point",default="", type=str,help="Mass point to process")
parser.add_argument("--local", dest="local",default=False, action='store_true', help="Process local files")

args, unknown = parser.parse_known_args()
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()

figdir = "./figures/signal_unc_%.0f_%.0f/" % (float(args.xgb_min)*100., float(args.xgb_max)*100.)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

if args.local:
    path="./trees/"
    sgn_masspoints=["500"]
    sgn_masspoint_files={
        "500" :"forMeas_bdt_v7_emu_scan_Zprime_M500_mcRun2018.root",
    }
else:
    path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
    sgn_masspoints=["200","400","600","800","1000"]
    sgn_masspoint_files={
        "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
        "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
        "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
        "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
        "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}

# Define the uncertainties
if args.local:
    uncs = [] # Not saved in the local processing
    scales = ['MuonES', 'ElectronES', 'JER', 'JES']
else:
    uncs = ['SFbtag', 'Electron_dxydz', 'Muon_dxydz', 'Trg', 'Prefire', 'PU', 'Electron_IsoID', 'Muon_IsoID',
            'Electron_ID', 'Muon_ID', 'Electron_RecoID', 'Muon_RecoID',
            'PtZ', 'PDF', 'Scale']
    scales = [] # Not implemented here


# Selection cuts
cuts = "(Flag_met && Flag_muon && Flag_electron) && xgb > %s && xgb <= %s" % (args.xgb_min, args.xgb_max)

# Loop through the signal files
effects = []
for mpoint in sgn_masspoints:
    if args.mass_point != "" and args.mass_point != mpoint: continue
    print "Processing mass point", mpoint
    # Create a TChain for the signal
    cc=rt.TChain("mytreefit")
    cc.Add(path+sgn_masspoint_files[mpoint])
    # Loop through each uncertainty
    for unc in uncs:
        [nom, up, down] = eval_unc(cc, unc, cuts, int(mpoint))
        effects.append([up/nom, down/nom])

    # Loop through the scale uncertainties
    for unc in scales:
        tree_up=rt.TChain("mytreefit")
        up_name = sgn_masspoint_files[mpoint]
        up_name = up_name.replace('mcRun', unc+'_up_mcRun')
        tree_up.Add(path+up_name)
        tree_down=rt.TChain("mytreefit")
        down_name = sgn_masspoint_files[mpoint]
        down_name = down_name.replace('mcRun', unc+'_down_mcRun')
        tree_down.Add(path+down_name)
        [nom, up, down] = eval_scale_unc(cc, tree_up, tree_down, unc, cuts, float(mpoint), figdir)
        effects.append([up/nom, down/nom])        
        
    for index in range(len(uncs)):
        print "%-15s: %.3f/%.3f" % (uncs[index], effects[index][0], effects[index][1])
    for index in range(len(scales)):
        print "%-15s: %.3f/%.3f" % (scales[index], effects[index+len(uncs)][0], effects[index+len(uncs)][1])
