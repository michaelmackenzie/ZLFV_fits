# Process Z->e+mu group uncertainties
import os
import argparse
import ROOT as rt
from array import array
from math import sqrt
from multiprocessing import Process

def extract_results(name):
    file_name = 'fitDiagnostics.' + name + '.root'
    f = rt.TFile.Open(file_name)
    if not f: return [-1., 0., 0.]
    t = f.Get('tree_fit_sb')
    t.GetEntry(0)
    return [t.r, t.rLoErr, t.rHiErr]

    
parser = argparse.ArgumentParser()
parser.add_argument("--card", dest="card",default="datacard_prm_units_v5_total.txt", type=str,help="datacard name")
parser.add_argument("--tag", dest="tag",default="test", type=str,help="output tag")
parser.add_argument("--rrange", dest="rrange",default="5", type=str,help="POI r range")
parser.add_argument("--skip-fits", dest="skip_fits",action='store_true', default=False,help="Skip fits")
args, unknown = parser.parse_known_args()
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit(1)


groups = {
    "envelope" : ["pdfindex_bin1", "pdfindex_bin2", "pdfindex_bin3"],
    "params" : [],
    "IDs_trig" : ["idEleSF", "idMuSF", "isoEleSF", "isoMuSF", "dxyEleSF", "dxyMuSF", "recoEleSF", "triggersSF"],
    "mu_scale" : ["mu_scale_bin1", "mu_scale_bin2", "mu_scale_bin3"],
    "ele_scale" : ["ele_scale_bin1", "ele_scale_bin2", "ele_scale_bin3"],
    "stats" : ["statsBin1", "statsBin2", "statsBin3"],
    "zmumu_yield" : ["Zmm"],
    "jer_jes" : ["jer","jes"],
    "scales" : ["bdt", "btagSF", "puSF", "prefire", "mixSF", "xsec", "zptSF", "pdf"],
    "all" : [],
}

results = {}

base_dir = os.getenv('CMSSW_BASE') + '/src/ZLFV_fits'

if not args.skip_fits:
    print ">>> Processing the nominal fit"
    command = base_dir + '/tests/zemu_fit.sh ' + args.card + ' "-t -1" ' + args.tag + '_nominal ' + args.rrange
    os.system(command)
results["nominal"] = extract_results(args.tag + '_nominal')

for group in groups:
    params = groups[group]
    if not args.skip_fits:
        command = base_dir + '/tests/zemu_fit.sh ' + args.card + ' "-t -1 --freezeParameters '
        if group == "params":
            command += 'var{bkg_.*},var{ratio_bkg_.*},rgx{pdfindex_.*}'
        elif group == "all":
            command += 'var{bkg_.*},var{ratio_bkg_.*},rgx{.*}'
        else:
            for param in params: command += param + ','
        command += '" ' + args.tag + '_' + group + ' ' + args.rrange
        print ">>> Processing group fit for", group
        os.system(command)
    results[group] = extract_results(args.tag + '_' + group)

print "Systematics omitted :     r     up    down  eff_up eff_down (rel_up, rel_down)"

n_res = results["nominal"]
print "%-20s: %7.3f +%.3f -%.3f" % ("nominal", n_res[0], n_res[2], n_res[1])
for group in groups:
    g_res = results[group]
    effects = [n_res[0] - g_res[0], sqrt(max(0., n_res[1]**2 - g_res[1]**2)), sqrt(max(0., n_res[2]**2 - g_res[2]**2))]
    print "%-20s: %7.3f +%.3f -%.3f +%.3f -%.3f (+%.3f -%.3f)" % (group, g_res[0], g_res[2], g_res[1], effects[2], effects[1], effects[2]/n_res[2], effects[1]/n_res[1])
    if group == "all":
        print "%-20s: %7.3f +%.3f -%.3f               (+%.3f -%.3f)" % ("Statistical", g_res[0], g_res[2], g_res[1], g_res[2]/n_res[2], g_res[1]/n_res[1])
