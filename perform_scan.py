# Process COMBINE datacards for Z'->emu scan
import os
import argparse
import ROOT as rt
from array import array

#----------------------------------------------------------------------------------------
# Process a single mass point
def process_datacard(card, directory, name, asimov = False, skip_fit = False, verbose = 0):
   if not os.path.isfile(directory + card):
      print "Card %s not found" % (card)
   if verbose > -1: print 'Processing mass point', name, ' (card =', card, ')'

   #----------------------------------------------------------------------------
   # Perform a signal rate fit
   #----------------------------------------------------------------------------

   if not skip_fit:
      command = 'combine -d %s -n .%s -M FitDiagnostics' % (card, name)
      if asimov: command += ' -t -1'
      #Allow negative measured signal rates
      command += ' --rMin -200 --rMax 200'
      #Additional commands to help fit converge properly
      command += ' --cminDefaultMinimizerStrategy 0'
      command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
      command += ' --cminRunAllDiscreteCombinations'
      output = 'fit_rate_%s.log' % (name)
      if verbose > 1: print command
      os.system('cd %s; %s >| %s; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Perform an upper limit evaluation
   #----------------------------------------------------------------------------

   if not skip_fit:
      command = 'combine -d %s -n .%s -M AsymptoticLimits' % (card, name)
      if asimov: command += ' -t -1'
      #Allow negative measured signal rates
      command += ' --rMin -200 --rMax 200'
      #Additional commands to help fit converge properly
      command += ' --cminDefaultMinimizerStrategy 0'
      command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
      command += ' --cminRunAllDiscreteCombinations'
      output = 'fit_limit_%s.log' % (name)
      if verbose > 1: print command
      os.system('cd %s; %s >| %s; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Extract the mass FIXME: Make this more robust
   #----------------------------------------------------------------------------

   ws_name = name
   ws_name = name.replace('bdt_', 'bdt_0d7_1d0_')
   ws_file = 'WorkspaceScanSGN/workspace_scansgn_v2_%s.root' % (ws_name)
   f = rt.TFile.Open(ws_file, 'READ')
   ws = f.Get('workspace_signal')
   mass_var = ws.var('mean_bin2')
   mass = mass_var.getVal()
   f.Close()
   

   #----------------------------------------------------------------------------
   # Extract the results
   #----------------------------------------------------------------------------

   fit_file = '%sfitDiagnostics.%s.root' % (directory, name)
   f = rt.TFile.Open(fit_file, 'READ')
   t = f.Get('tree_fit_sb')
   t.GetEntry(0)
   r_fit = [t.r, t.rLoErr, t.rHiErr]
   if verbose > 0: print 'r fit results:', r_fit
   f.Close()

   lim_file = '%shiggsCombine.%s.AsymptoticLimits.mH120.root' % (directory, name)
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

   return [r_fit, r_lim, mass]


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="bdt_v01", type=str,help="datacard directory name")
parser.add_argument("--skip-fit", dest="skip_fit",default=False, action='store_true',help="Skip fits, assume already processed")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="Perform fits to the data instead of Asimov fits")

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
figdir = "./figures/scan_%s/" % (args.name)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))
   
rt.gROOT.SetBatch(True)

#----------------------------------------------
# Perform the scan
#----------------------------------------------

list_of_files = os.listdir(carddir)
asimov = not args.unblind

# List of results
masses = array('d')
masses_errs = array('d')
r_fits = array('d')
r_lo_errs = array('d')
r_hi_errs = array('d')
r_lims = array('d')
r_exps = array('d')
r_exps_lo = array('d') #1 sigma band
r_exps_hi = array('d')
r_exps_lo_2 = array('d') #2 sigma band
r_exps_hi_2 = array('d')

prev_mass = -1.
min_lim =  1.e10
max_lim = -1.e10
min_r =  1.e10
max_r = -1.e10
for f in list_of_files:
   if '.txt' not in f: continue
   # only process the merged fits
   if '0d7' in f: continue
   mass_point = f.split('_mp')[1].split('.txt')[0]
   [r_fit, r_lim, mass] = process_datacard(f, carddir, args.name + '_mp'+mass_point, asimov, args.skip_fit)

   # store the results
   masses.append(mass)
   masses_errs.append((mass - prev_mass)/2. if prev_mass > 0. else 0.5) #rough mass bin estimate
   r_fits.append(r_fit[0])
   r_lo_errs.append(r_fit[1])
   r_hi_errs.append(r_fit[2])
   r_lims.append(r_lim[-1])
   r_exps.append(r_lim[2])  
   r_exps_lo.append(r_lim[2] - r_lim[1])
   r_exps_hi.append(r_lim[3] - r_lim[2])
   r_exps_lo_2.append(r_lim[2] - r_lim[0])
   r_exps_hi_2.append(r_lim[4]- r_lim[2])

   prev_mass = mass
   min_lim = min(r_lim[0], r_lim[-1], min_lim)
   max_lim = max(r_lim[4], r_lim[-1], max_lim)
   min_r = min(r_fit[0] - r_fit[1], min_r)
   max_r = max(r_fit[0] + r_fit[2], max_r)
 
#----------------------------------------------
# Plot the results
#----------------------------------------------

#----------------------------------------------
# Limit plot
#----------------------------------------------

g_exp   = rt.TGraph(len(masses), masses, r_exps)
g_exp_1 = rt.TGraphAsymmErrors(len(masses), masses, r_exps, masses_errs, masses_errs, r_exps_lo  , r_exps_hi  )
g_exp_2 = rt.TGraphAsymmErrors(len(masses), masses, r_exps, masses_errs, masses_errs, r_exps_lo_2, r_exps_hi_2)

c = rt.TCanvas('c_lim', 'c_lim', 800, 600)
g_exp_2.SetTitle("95% CL_{S} limit vs. Z' mass; Z' mass (GeV/c^{2}); #sigma(Z')*BR(Z'->e#mu) (fb)")
g_exp_2.SetFillColor(rt.kOrange)
g_exp_2.SetLineColor(rt.kOrange)
g_exp_1.SetFillColor(rt.kGreen+1)
g_exp_1.SetLineColor(rt.kGreen+1)
g_exp_2.Draw("AE3")
g_exp_1.Draw("E3")
g_exp.SetLineStyle(rt.kDashed)
g_exp.SetLineWidth(2)
g_exp.SetLineColor(rt.kBlack)
g_exp.Draw("XL")

g_obs = rt.TGraph(len(masses), masses, r_lims)
g_obs.SetMarkerStyle(20)
g_obs.SetMarkerSize(0.8)
g_obs.SetLineWidth(2)
g_obs.SetMarkerColor(rt.kBlack)
g_obs.SetLineColor(rt.kBlack)
g_obs.Draw("PL")

g_exp_2.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_exp_2.GetYaxis().SetRangeUser(0.8*min_lim, 1.4*max_lim)

leg = rt.TLegend(0.7, 0.7, 0.89, 0.89)
leg.SetLineWidth(0)
leg.AddEntry(g_obs  , 'Observed', 'PL')
leg.AddEntry(g_exp  , 'Expected', 'L')
leg.AddEntry(g_exp_1, '#pm1#sigma' , 'F')
leg.AddEntry(g_exp_2, '#pm2#sigma' , 'F')
leg.Draw()


c.SaveAs(figdir+'limits.png')

#----------------------------------------------
# Signal rate plot
#----------------------------------------------

g_r = rt.TGraphAsymmErrors(len(masses), masses, r_fits, masses_errs, masses_errs, r_lo_errs  , r_hi_errs)

c = rt.TCanvas('c_fit', 'c_fit', 800, 600)
g_r.SetTitle("Z' fit rate vs. Z' mass; Z' mass (GeV/c^{2}); #sigma(Z')*BR(Z'->e#mu) (fb)")
g_r.SetMarkerStyle(20)
g_r.SetMarkerSize(0.8)
g_r.SetLineWidth(2)
g_r.SetMarkerColor(rt.kBlack)
g_r.SetLineColor(rt.kBlack)
g_r.Draw("APE1")

g_r.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_r.GetYaxis().SetRangeUser(min_r - 0.05*(max_r - min_r), max_r + 0.1*(max_r - min_r))
c.SaveAs(figdir+'fits.png')
