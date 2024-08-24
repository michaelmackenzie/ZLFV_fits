# Compare two identical mass point COMBINE Z'->emu scan sensitivities
import os
import argparse
import ROOT as rt
from array import array
from math import exp
from math import log10

#----------------------------------------------------------------------------------------
# Define the sorting of the datacards
def file_sort(f):
   return int(f.split('_mp')[1].replace('.txt',''))

#----------------------------------------------------------------------------------------
# Process a single mass point
def process_datacard(card, directory, name, tag, verbose = 0):
   if verbose > -1: print 'Processing mass point', name, '(tag =', tag+')'

   #----------------------------------------------------------------------------
   # Extract the mass FIXME: Make this more robust
   #----------------------------------------------------------------------------

   ws_name = name
   ws_name = name.replace('bdt_', 'bdt_0d7_1d0_')
   ws_file = 'WorkspaceScanSGN/workspace_scansgn_v2_%s.root' % (ws_name)
   if 'toy' in ws_file:
      ws_file = ws_file.split('_toy')[0] + '_mp' + ws_file.split('_mp')[1]
   f = rt.TFile.Open(ws_file, 'READ')
   ws = f.Get('workspace_signal')
   mass_var = ws.var('mean_bin2')
   mass = mass_var.getVal()
   f.Close()
   

   #----------------------------------------------------------------------------
   # Extract the results
   #----------------------------------------------------------------------------
   
   lim_file = '%shiggsCombine.%s%s_asimov.AsymptoticLimits.mH120.root' % (directory, name, tag)
   f = rt.TFile.Open(lim_file, 'READ')
   t = f.Get('limit')
   r_lim = []
   for index in range(t.GetEntries()):
      t.GetEntry(index)
      r_lim.append(t.limit)
   if verbose > 0: print 'r limit results:', r_lim
   f.Close()
   
   #----------------------------------------------------------------------------
   # Return the results
   #----------------------------------------------------------------------------

   return [r_lim, mass]


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="bdt_v01", type=str,help="datacard directory name")
parser.add_argument("--max-steps", dest="max_steps",default=-1, type=int, help="Maximum steps to take in the scan")
parser.add_argument("--first-step", dest="first_step",default=0, type=int, help="First mass step to process")
parser.add_argument("--card-tag-1", dest="card_tag_1",type=str, help="Card name tag for version 1", required=True)
parser.add_argument("--card-tag-2", dest="card_tag_2",type=str, help="Card name tag for version 2", required=True)
parser.add_argument("--tag", dest="tag",default="", type=str, help="Output directory tag")
parser.add_argument("-v", dest="verbose",default=0, type=int,help="Add verbose printout")

args, unknown = parser.parse_known_args()


##### configuaration ########

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()

#----------------------------------------------
# Setup the scan
#----------------------------------------------

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/scan_ratio_%s%s/" % (args.name, args.tag)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
   
rt.gROOT.SetBatch(True)

#----------------------------------------------
# Perform the scan
#----------------------------------------------

list_of_files = [f for f in os.listdir(carddir) if '.txt' in f and '0d7' not in f]


# Sort the list by mass point
list_of_files.sort(key=file_sort)
if args.first_step > 0:
   list_of_files = list_of_files[args.first_step:]
if args.max_steps > 0:
   list_of_files = list_of_files[:args.max_steps]


# List of results
masses = array('d')
r_exp_ratios = array('d')


min_ratio =  1.e10
max_ratio = -1.e10
for f in list_of_files:
   mass_point = f.split('_mp')[1].split('.txt')[0]
   [r_lim_1, mass_1] = process_datacard(f, carddir, args.name + '_mp'+mass_point, args.card_tag_1, args.verbose)
   [r_lim_2, mass_2] = process_datacard(f, carddir, args.name + '_mp'+mass_point, args.card_tag_2, args.verbose)
   if abs(mass_1-mass_2) > 1.e-4:
      print "Masses not equal! Mass 1 = %.4f, mass 2 = %.4f" % (mass_1, mass_2)

   # store the results
   masses.append(mass_1)
   r_exp_ratios.append(r_lim_2[2] / r_lim_1[2])

   min_ratio = min(r_exp_ratios[-1], min_ratio)
   max_ratio = max(r_exp_ratios[-1], max_ratio)


#----------------------------------------------
# Plot the results
#----------------------------------------------

rt.gStyle.SetOptStat(0)

#----------------------------------------------
# Limit plot
#----------------------------------------------

c = rt.TCanvas('c_lim', 'c_lim', 800, 600)
g = rt.TGraph(len(masses), masses, r_exp_ratios)
g.SetTitle("95% CL_{S} limit ratio vs. Z' mass; Z' mass (GeV/c^{2}); #sigma(Z')*BR(Z'->e#mu) ratio")
g.SetMarkerStyle(20)
g.SetMarkerSize(0.8)
g.SetLineWidth(2)
g.SetMarkerColor(rt.kBlack)
g.SetLineColor(rt.kBlack)
g.Draw("APL")

g.GetXaxis().SetRangeUser(masses[0], masses[-1])
g.GetYaxis().SetRangeUser(0.8*min_ratio, 1.2*max_ratio)
c.SaveAs(figdir+'limit_ratios.png')
