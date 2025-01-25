# Read in toy scan results and create a distribution of minimum p-values

import os
import argparse
import ROOT as rt
from array import array

def retrieve_info(card, directory, name, tag = '', verbose = 0):

   #----------------------------------------------------------------------------
   # Extract the results
   #----------------------------------------------------------------------------

   mass = float(card.split('mass-')[1].split('_')[0])

   fit_file = '%sfitDiagnostics.%s%s.root' % (directory, name, tag)
   f = rt.TFile.Open(fit_file, 'READ')
   if not f: return [None, None, None, mass]
   t = f.Get('tree_fit_sb')
   try:
       t.GetEntry(0)
   except:
       print "Tree not found in fit file", fit_file
       return [None, None, None, mass]
   r_fit = [t.r, t.rLoErr, t.rHiErr]
   if verbose > 0: print 'r fit results:', r_fit
   f.Close()

   sig_file = '%shiggsCombine.%s%s.Significance.mH120.root' % (directory, name, tag)
   f = rt.TFile.Open(sig_file, 'READ')
   if not f: return [None, None, None, mass]
   t = f.Get('limit')
   try:
       t.GetEntry(0)
   except:
       print "Tree not found in sig file", sig_file
       return [None, None, None, mass]
   r_sig = t.limit
   if verbose > 0: print 'r significance results:', r_sig
   f.Close()
   
   lim_file = '%shiggsCombine.%s%s.AsymptoticLimits.mH120.root' % (directory, name, tag)
   f = rt.TFile.Open(lim_file, 'READ')
   if not f: return [None, None, None, mass]
   t = f.Get('limit')
   r_lim = []
   try:
       for index in range(t.GetEntries()):
           t.GetEntry(index)
           r_lim.append(t.limit)
   except:
       print "Tree entries not found in limit file", lim_file
       return [None, None, None, mass]
   if verbose > 0: print 'r limit results:', r_lim
   f.Close()
   
   #----------------------------------------------------------------------------
   # Return the results
   #----------------------------------------------------------------------------

   return [r_fit, r_lim, r_sig, mass]


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="bdt_v01", type=str,help="datacard directory name")
parser.add_argument("--ntoys", dest="ntoys",default=1, type=int, help="Number of toys processed")
parser.add_argument("--card-tag", dest="card_tag",default="", type=str, help="Card name tag to process")
parser.add_argument("--tag", dest="tag",default="", type=str, help="Output directory tag")
parser.add_argument("-v", '--verbose', dest="verbose",default=0, type=int,help="Add verbose printout")

args, unknown = parser.parse_known_args()

# For each toy, get the list of measurement results, store the most extreme p-value
min_pval = array('f', [0])
max_sig  = array('f', [0])
min_sig  = array('f', [0])
os.system('[ ! -d lee ] && mkdir lee')
fout = rt.TFile('lee/lee_effect_%s.root' % (args.name), 'RECREATE')
tree = rt.TTree('LEE', 'Look-elsewhere effect Tree')
tree.Branch('min_pvalue', min_pval, 'min_pvalue/F')
tree.Branch('max_sig'   , max_sig , 'max_sig/F'   )
tree.Branch('min_sig'   , min_sig , 'min_sig/F'   )
for toy in range(1,args.ntoys+1):
    # Get a list of the datacards
    carddir = "./datacards/%s_toy_%i/" % (args.name, toy)
    list_of_files = [f for f in os.listdir(carddir) if '.txt' in f]
    # If not using a BDT score region tag, only take the merged files
    if args.card_tag == "":  list_of_files = [f for f in list_of_files if '0d7' not in f]
    else:                    list_of_files = [f for f in list_of_files if args.card_tag in f]

    if len(list_of_files) == 0:
        print "No card files found for toy %i!" % (toy)
        break

    if len(list_of_files) < 70:
        print "Toy %2i is likely missing cards, only has %i cards found" % (toy, len(list_of_files))
        continue
    min_pval[0] = 1.
    max_sig [0] = -10.
    min_sig [0] =  10.
    nfail = 0
    for f in list_of_files:
        mass_point = f.split('_mp')[1].split('.txt')[0]
        [r_fit, r_lim, r_sig, mass] = retrieve_info(f, carddir, args.name+'_toy_'+str(toy)+'_mp'+mass_point, args.tag, max(0, args.verbose-2))
        if r_fit is None or r_lim is None or r_sig is None:
            nfail+=1
            continue
        pval = rt.RooStats.SignificanceToPValue(r_sig)
        min_pval[0] = min(pval , min_pval[0])
        max_sig [0] = max(r_sig, max_sig [0])
        min_sig [0] = min(r_sig, min_sig [0])
        if args.verbose > 1: print " Mass hypothesis %5.1f: p = %.4f, sig = %.2f" % (mass, pval, r_sig)
    if nfail > 0.1*len(list_of_files):
        print "Toy %2i had %i failures, skipping this toy" % (toy, nfail)
    if args.verbose > 0: print "Toy %3i: min_pval = %.5f, max_sig = %.2f, min_sig = %.2f" % (toy, min_pval[0], max_sig[0], min_sig[0])
    tree.Fill()


# Make a plot of each value
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)
figdir = "./figures/scan_lee_%s%s/" % (args.name, args.tag)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

c = rt.TCanvas()
h = rt.TH1F('h_pval', 'LEE log10(p-value) distribution', 100, -6., 0.)
tree.Draw("log10(min_pvalue) >> h_pval")
h.SetLineWidth(2)
h.GetXaxis().SetTitle('log10(minimum p)')
c.SaveAs(figdir + 'min_pvals.png')

h = rt.TH1F('h_sig', 'LEE minimum significance distribution', 40, -6., 0.)
tree.Draw("min_sig >> h_sig")
h.SetLineWidth(2)
h.GetXaxis().SetTitle('minimum significance')
c.SaveAs(figdir + 'min_sigs.png')

h = rt.TH1F('h_max_sig', 'LEE maximum significance distribution', 40, 0., 6.)
tree.Draw("max_sig >> h_max_sig")
h.SetLineWidth(2)
h.GetXaxis().SetTitle('maximum significance')
c.SaveAs(figdir + 'max_sigs.png')

fout.cd()
tree.Write()
fout.Close()
