# Process COMBINE datacards for the Z'->emu scan
import os
import argparse
import ROOT as rt
from array import array
from math import exp
from math import log10
from math import sqrt
from multiprocessing import Process

#----------------------------------------------------------------------------------------
# Define the sorting of the datacards
def file_sort(f):
   return int(f.split('_mp')[1].replace('.txt',''))

#----------------------------------------------------------------------------------------
# Draw the CMS label
def draw_cms_label(is_prelim = True):
    #CMS prelim drawing
    cmslabel = rt.TText()
    cmslabel.SetNDC();
    cmslabel.SetTextColor(1);
    cmslabel.SetTextSize(0.06);
    cmslabel.SetTextAlign(11);
    cmslabel.SetTextAngle(0);
    cmslabel.SetTextFont(61);
    cmslabel.DrawText(0.14, 0.83, "CMS");
    if is_prelim:
      cmslabel.SetTextFont(52);
      cmslabel.SetTextSize(0.76*cmslabel.GetTextSize());
      cmslabel.DrawText(0.14, 0.78, "Preliminary");

def draw_luminosity(year = -1):
  label = rt.TLatex()
  label.SetNDC()
  label.SetTextFont(42)
  label.SetTextSize(0.04)
  label.SetTextAlign(31)
  label.SetTextAngle(0)
  period = "%i, " % (year) if year > 2000 else ""
  lums = [36.33, 41.48, 59.83, 137.64]
  lum = lums[year-2016] if year > 2000 else lums[-1]
  label.DrawLatex(0.97, 0.915, "%s%.0f fb^{-1} (13 TeV)" % (period,lum));

#----------------------------------------------------------------------------------------
# Smooth expected limits by fitting to an exponential in a moving window
def alt_smooth_limits(vals, masses):
   func = rt.TF1('func', 'exp([0] + [1]*(x/1000))', masses[0], masses[-1])
   func.SetParameters(1., -1.)
   g = rt.TGraph(len(masses), masses, vals)
   # Fit the graph in a moving window
   for index in range(len(masses)):
      mass = masses[index]
      min_mass = masses[max(0, index-6)]
      max_mass = masses[min(len(masses)-1, index+6)]
      func.SetRange(min_mass, max_mass)
      g.Fit(func, 'R Q 0')
      vals[index] = func.Eval(mass)
   # Make sure it's monotomic if multiple fits are used
   for index in range(1, len(masses)-1):
      if ((vals[index] < vals[index-1] and vals[index] < vals[index+1]) or
          (vals[index] > vals[index-1] and vals[index] > vals[index+1])):
         vals[index] = 0.5*(vals[index+1] + vals[index-1])
      
   return vals

#----------------------------------------------------------------------------------------
# Smooth expected limits by fitting to an exponential
def smooth_limits(vals, masses):
   mode = 5 # How to smooth the fits
   if mode == 0:
      func = rt.TF1('func', 'exp([0] + [1]*(x/1000))', masses[0], masses[-1])
      func.SetParameters(1., -1.)
   elif mode == 1:
      func = rt.TF1('func', 'exp([0] + [1]*(x/1000)) + exp([2] + [3]*(x/1000))', masses[0], masses[-1])
      func.SetParameters(1., -1., 0., -0.9)
   elif mode == 2:
      func = rt.TF1('func', '[0] + [1]*(x/1000) + [2]*pow((x/1000),2) + [3]*pow((x/1000),3)', masses[0], masses[-1])
      func.SetParameters(10., -1., 0., 0.) 
   elif mode == 3:
      func = rt.TF1('func', 'pow([0]*(x/1000),[1]) + exp([2] + [3]*(x/1000))', masses[0], masses[-1])
      func.SetParameters(3., -1.3, 0., -0.9)
   elif mode == 4: # Fit slightly away from the start
      func = rt.TF1('func', 'exp([0] + [1]*(x/1000))', 120., masses[-1])
      func.SetParameters(1., -1.)
   elif mode == 5: # Fit in different regions
      func = rt.TF1('func', '(x < 120.) * exp([0] + [1]*(x/1000)) + (x > 120)*exp([2] + [3]*(x/1000))', masses[0], masses[-1])
      func.SetParameters(1., -1., 1., -1.)
      # func = rt.TF1('func', '(x < 120.) * exp([0] + [1]*(x/1000)) + (x > 120) * (x < 130)*exp([2] + [3]*(x/1000)) + (x > 130)*exp([4] + [5]*(x/1000))', masses[0], masses[-1])
      # func.SetParameters(1., -1., 1., -1., 1., -1.)
   else:
      print "Unknown smoothing mode", mode, "--> no smoothing will be performed!"
      return vals
   g = rt.TGraph(len(masses), masses, vals)
   g.Fit(func, 'R Q 0')
   for index in range(len(masses)):
      mass = masses[index]
      if mass < func.GetXmin(): continue
      if mass > func.GetXmax(): continue
      vals[index] = func.Eval(mass)
   if mode == 5:
      # Make sure it's monotomic if multiple fits are used
      for index in range(1, len(masses)-1):
         if ((vals[index] < vals[index-1] and vals[index] < vals[index+1]) or
             (vals[index] > vals[index-1] and vals[index] > vals[index+1])):
            vals[index] = 0.5*(vals[index+1] + vals[index-1])
      
   return vals

#----------------------------------------------------------------------------------------
# Process a single mass point
def process_datacard(card, directory, name, asimov = False, tag = '', verbose = 0):
   if not os.path.isfile(directory + card):
      print "Card %s not found" % (card)
   if verbose > -1: print 'Processing mass point', name, '(card =', card+')'

   #----------------------------------------------------------------------------
   # Perform a signal rate fit
   #----------------------------------------------------------------------------

   if asimov: tag += '_asimov'
   command = 'combine -d %s -n .%s%s -M FitDiagnostics' % (card, name, tag)
   if asimov: command += ' -t -1'
   #Allow negative measured signal rates
   command += ' --rMin -50 --rMax 50'
   #Additional commands to help fit converge properly
   command += ' --cminDefaultMinimizerStrategy 0'
   command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
   command += ' --cminRunAllDiscreteCombinations'
   # command += '  --X-rtd MINIMIZER_multiMin_maskChannels=2'
   command += ' --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1'
   command += ' --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001'
   # command += ' --rAbsAcc 0.0005 --rRelAcc 0.0005'
   output = 'fit_rate_%s.log' % (name)
   if verbose > 1: print command
   os.system('cd %s; %s >| %s 2>&1; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Evaluate the signal rate significance
   #----------------------------------------------------------------------------

   command = 'combine -d %s -n .%s%s -M Significance --uncapped 1' % (card, name, tag)
   if asimov: command += ' -t -1'
   #Allow negative measured signal rates
   command += ' --rMin -50 --rMax 50'
   #Additional commands to help fit converge properly
   command += ' --cminDefaultMinimizerStrategy 0'
   command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
   command += ' --cminRunAllDiscreteCombinations'
   # command += ' --X-rtd MINIMIZER_multiMin_maskChannels=2'
   command += ' --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1'
   command += ' --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001'
   # command += ' --rAbsAcc 0.0005 --rRelAcc 0.0005'
   output = 'fit_sig_%s.log' % (name)
   if verbose > 1: print command
   os.system('cd %s; %s >| %s 2>&1; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Perform an upper limit evaluation
   #----------------------------------------------------------------------------

   command = 'combine -d %s -n .%s%s -M AsymptoticLimits' % (card, name, tag)
   if asimov: command += ' -t -1 --noFitAsimov' # --freezeParameters pdfindex_bin1,pdfindex_bin2
   #Allow negative measured signal rates
   command += ' --rMin -50 --rMax 50'
   #Additional commands to help fit converge properly
   command += ' --cminDefaultMinimizerStrategy 0'
   command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
   command += ' --cminRunAllDiscreteCombinations'
   # command += ' --X-rtd MINIMIZER_multiMin_maskChannels=2'
   command += ' --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1'
   command += ' --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001'
   command += ' --rAbsAcc 0.0005 --rRelAcc 0.0005'
   output = 'fit_limit_%s.log' % (name)
   if verbose > 1: print command
   os.system('cd %s; %s >| %s 2>&1; cd ..' % (directory, command, output))

#----------------------------------------------------------------------------------------
# Retrieve fit information for a single mass point
def retrieve_info(card, directory, name, asimov = False, tag = '', verbose = 0):

   mass = float(card.split('mass-')[1].split('_')[0])
   if asimov: tag += '_asimov'

   #----------------------------------------------------------------------------
   # Extract the results
   #----------------------------------------------------------------------------

   fit_file = '%sfitDiagnostics.%s%s.root' % (directory, name, tag)
   f = rt.TFile.Open(fit_file, 'READ')
   t = f.Get('tree_fit_sb')
   t.GetEntry(0)
   r_fit = [t.r, t.rLoErr, t.rHiErr]
   if verbose > 0: print 'r fit results:', r_fit
   f.Close()

   sig_file = '%shiggsCombine.%s%s.Significance.mH120.root' % (directory, name, tag)
   f = rt.TFile.Open(sig_file, 'READ')
   t = f.Get('limit')
   t.GetEntry(0)
   r_sig = t.limit
   if verbose > 0: print 'r significance results:', r_sig
   f.Close()
   
   lim_file = '%shiggsCombine.%s%s.AsymptoticLimits.mH120.root' % (directory, name, tag)
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

   return [r_fit, r_lim, r_sig, mass]


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="bdt_v01", type=str,help="datacard directory name")
parser.add_argument("-j", "--nthreads", dest="nthreads",default=8,type=int,help="Number of threads to process using")
parser.add_argument("--skip-fits", dest="skip_fits",default=False, action='store_true',help="Skip fits, assume already processed")
parser.add_argument("--asimov", dest="asimov",default=False, action='store_true',help="Perform fits Asimov dataset")
parser.add_argument("--smooth-expected", dest="smooth_expected",default=False, action='store_true',help="Smooth the expected limit distribution")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="Plot the observed limits")
parser.add_argument("--draw-asimov", dest="draw_asimov",default=False, action='store_true',help="Draw the observed fits in the Asimov scan")
parser.add_argument("--max-steps", dest="max_steps",default=-1, type=int, help="Maximum steps to take in the scan")
parser.add_argument("--first-step", dest="first_step",default=0, type=int, help="First mass step to process")
parser.add_argument("--card-tag", dest="card_tag",default="", type=str, help="Card name tag to process")
parser.add_argument("--tag", dest="tag",default="", type=str, help="Output directory tag")
parser.add_argument("--extra-tag", dest="extra_tag",default="", type=str, help="Extra tag on output figure directory")
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

if args.tag != "": args.tag = "_" + args.tag
if args.extra_tag != "": args.extra_tag = "_" + args.extra_tag

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/scan_%s%s%s%s/" % (args.name, "_asimov" if args.asimov else "", args.tag, args.extra_tag)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))
   
rt.gROOT.SetBatch(True)

# Get the local --> global p-value TTree if available
lee_file = rt.TFile.Open("lee/lee_effect_%s.root" % (args.name), "READ")
lee_tree = lee_file.Get('LEE') if lee_file else None
n_lee_tree = lee_tree.GetEntries() if lee_tree else 0

#----------------------------------------------
# Perform the scan
#----------------------------------------------

list_of_files = [f for f in os.listdir(carddir) if '.txt' in f]
# If not using a BDT score region tag, only take the merged files
if args.card_tag == "":  list_of_files = [f for f in list_of_files if '0d7' not in f]
else:                    list_of_files = [f for f in list_of_files if args.card_tag in f]

if len(list_of_files) == 0:
   print "No card files found!"
   exit()

# Sort the list by mass point
list_of_files.sort(key=file_sort)
if args.first_step > 0:
   list_of_files = list_of_files[args.first_step:]
if args.max_steps > 0:
   list_of_files = list_of_files[:args.max_steps]

asimov = args.asimov

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
r_sigs = array('d')
global_pvals = array('d')


prev_mass = -1.
min_lim =  1.e10
max_lim = -1.e10
min_r =  1.e10
max_r = -1.e10

jobs = [] # For multithreaded processing

# Process the combine jobs for each card
if not args.skip_fits:
   for f in list_of_files:
      if '.txt' not in f: continue
      mass_point = f.split('_mp')[1].split('.txt')[0]
      job = Process(target = process_datacard, args=(f, carddir, args.name+ '_mp'+mass_point, asimov, args.tag, args.verbose))
      jobs.append(job)
   if args.nthreads < 1: args.nthreads = 1
   print("Parallel processing using %i threads" % (args.nthreads))
   for ithread in range(0,len(jobs),args.nthreads):
      nthread = args.nthreads+ithread
      if (nthread>len(jobs)):
         nthread=len(jobs)
      print("Processing threads %i to %i out of %i jobs: %5.1f%% processed" % (ithread, nthread-1, len(jobs), (ithread*100./len(jobs))))
      for job in jobs[ithread:nthread]:
         job.start()
      for job in jobs[ithread:nthread]:    
         job.join()

# Retrieve and plot the fit results for each card
for f in list_of_files:
   if '.txt' not in f: continue
   mass_point = f.split('_mp')[1].split('.txt')[0]
   [r_fit, r_lim, r_sig, mass] = retrieve_info(f, carddir, args.name+ '_mp'+mass_point, asimov, args.tag, args.verbose)

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
   r_sigs.append(r_sig)

   prev_mass = mass
   min_lim = min(r_lim[0], r_lim[-1], min_lim)
   max_lim = max(r_lim[4], r_lim[-1], max_lim)
   min_r = min(r_fit[0] - r_fit[1], min_r)
   max_r = max(r_fit[0] + r_fit[2], max_r)

   # Print the LEE information
   if lee_tree and n_lee_tree > 0:
      n_more_extreme = lee_tree.Draw("", "max_sig > %.3f" % (r_sig), "goff")
      global_p = n_more_extreme * 1. / n_lee_tree
      global_unc = sqrt(global_p*(1.-global_p)/n_lee_tree)
      local_p = rt.RooStats.SignificanceToPValue(r_sig)
      print "mass = %.1f, r = %4.1f, sigma = %5.2f, local p = %.4f, global p = %.4f +- %.4f%s" % (mass, r_fit[0], r_sig, local_p,
                                                                                                  global_p, global_unc,
                                                                                                  " <--" if r_sig > 2. else "")
      global_pvals.append(global_p)

if len(masses) > 1: masses_errs[0] = (masses[1] - masses[0])/2.

#----------------------------------------------
# Smooth the expected limits if requested
#----------------------------------------------

if args.smooth_expected:
   # Fit the expected limits/uncertainties with an exponential
   r_exps      = alt_smooth_limits(r_exps     , masses)
   r_exps_lo   = alt_smooth_limits(r_exps_lo  , masses)
   r_exps_hi   = alt_smooth_limits(r_exps_hi  , masses)
   r_exps_lo_2 = alt_smooth_limits(r_exps_lo_2, masses)
   r_exps_hi_2 = alt_smooth_limits(r_exps_hi_2, masses)
   

#----------------------------------------------
# Plot the results
#----------------------------------------------

rt.gStyle.SetOptStat(0)
rt.gStyle.SetPadTickX(1);
rt.gStyle.SetPadTickY(1);

#----------------------------------------------
# Limit plot
#----------------------------------------------

# Limit plot colors
LimitGreen  = rt.TColor.GetColor("#607641")
LimitYellow = rt.TColor.GetColor("#F5BB54")

draw_obs = (args.asimov and args.draw_asimov) or args.unblind

g_exp   = rt.TGraph(len(masses), masses, r_exps)
g_exp_1 = rt.TGraphAsymmErrors(len(masses), masses, r_exps, masses_errs, masses_errs, r_exps_lo  , r_exps_hi  )
g_exp_2 = rt.TGraphAsymmErrors(len(masses), masses, r_exps, masses_errs, masses_errs, r_exps_lo_2, r_exps_hi_2)

c = rt.TCanvas('c_lim', 'c_lim', 800, 600)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.10)
c.SetBottomMargin(0.13)
g_exp_2.SetTitle('') #"95% CL_{S} limit vs. Z' mass"
g_exp_2.SetFillColor(LimitYellow)
g_exp_2.SetLineColor(LimitYellow)
g_exp_1.SetFillColor(LimitGreen)
g_exp_1.SetLineColor(LimitGreen)
g_exp_2.Draw("AE3")
g_exp_1.Draw("E3")
g_exp.SetLineStyle(rt.kDashed)
g_exp.SetLineWidth(2)
g_exp.SetLineColor(rt.kBlack)
g_exp.Draw("XL")
g_exp_2.GetXaxis().SetTitle("m_{Z\'} [GeV]");
g_exp_2.GetYaxis().SetTitle("#sigma(pp#rightarrowZ\'+X)#it{B}(Z\'#rightarrowe#mu) [fb]");

if draw_obs:
   g_obs = rt.TGraph(len(masses), masses, r_lims)
   g_obs.SetMarkerStyle(20)
   g_obs.SetMarkerSize(0.8)
   g_obs.SetLineWidth(2)
   g_obs.SetMarkerColor(rt.kBlack)
   g_obs.SetLineColor(rt.kBlack)
   g_obs.Draw("PL")

xmin = min(100., masses[0])
xmax = max(500., masses[-1])
ymin = 1.e-4
ymax = 1.1*max_lim if args.first_step > 0 else 1.03*max_lim
g_exp_2.GetXaxis().SetRangeUser(xmin, xmax)
g_exp_2.GetYaxis().SetRangeUser(ymin, ymax)
g_exp_2.GetXaxis().SetLabelSize(0.045)
g_exp_2.GetYaxis().SetLabelSize(0.045)
g_exp_2.GetXaxis().SetTitleSize(0.05)
g_exp_2.GetXaxis().SetTitleOffset(0.9)
g_exp_2.GetYaxis().SetTitleSize(0.05)
g_exp_2.GetYaxis().SetTitleOffset(0.6)

leg = rt.TLegend(0.6, 0.60, 0.89, 0.80)
leg.SetLineWidth(0)
if draw_obs:
   leg.AddEntry(g_obs  , 'Observed', 'PL')
leg.AddEntry(g_exp  , 'Expected', 'L')
leg.AddEntry(g_exp_1, '#pm1 std. deviation' , 'F')
leg.AddEntry(g_exp_2, '#pm2 std. deviation' , 'F')
leg.Draw()
leg.SetTextSize(0.04)

# Add limit type label
label = rt.TLatex()
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.04)
label.SetTextAlign(11)
label.SetTextAngle(0)
label.DrawLatex(0.6, 0.81, "95% CL upper limit")

is_prelim = False
draw_cms_label(is_prelim)
draw_luminosity()

c.SaveAs(figdir+'limits.pdf')
c.SaveAs(figdir+'limits.root')
g_exp_2.GetYaxis().SetRangeUser(0.5*min_lim, 5*max_lim)
c.SetLogy()
c.SaveAs(figdir+'limits_log.pdf')
c.SaveAs(figdir+'limits_log.root')

#----------------------------------------------
# Limit plot with log plot on top
#----------------------------------------------

c = rt.TCanvas('c_lim_2', 'c_lim_2', 800, 600)
pad1 = rt.TPad('pad1', 'pad1', 0., 0., 1., 1.)
# pad2 = rt.TPad('pad2', 'pad2', 0.3, 0.42, 0.7, 0.88)
pad2 = rt.TPad('pad2', 'pad2', 0.45, 0.30, 0.895, 0.88)
leg.SetX1NDC(0.2); leg.SetX2NDC(0.40);
pad2.SetTopMargin(0.01)
pad2.SetRightMargin(0.03)
pad2.SetBottomMargin(0.08)
pad2.SetLeftMargin(0.10)
pad1.Draw()
pad2.Draw()
pad1.cd()
g_exp_2.Draw("AE3")
g_exp_1.Draw("E3")
g_exp.Draw("XL")
if draw_obs:
   g_obs.Draw("PL")
leg.Draw()
g_exp_2.GetYaxis().SetRangeUser(0.*min_lim, 1.01*max_lim)
pad2.cd()
g_exp_2_log = g_exp_2.Clone('g_exp_2_log')
g_exp_2_log.Draw("AE3")
g_exp_1.Draw("E3")
g_exp.Draw("XL")
if draw_obs:
   g_obs.Draw("PL")
g_exp_2_log.SetTitle('')
# g_exp_2_log.GetXaxis().SetTitle("");
# g_exp_2_log.GetYaxis().SetTitle("");
g_exp_2_log.GetYaxis().SetRangeUser(0.5*min_lim, 1.7*max_lim)
pad2.SetLogy()
c.SaveAs(figdir+'limits_plus_log.png')



#----------------------------------------------
# Signal rate plot
#----------------------------------------------

if not args.asimov and not args.unblind: exit()

g_r = rt.TGraphAsymmErrors(len(masses), masses, r_fits, masses_errs, masses_errs, r_lo_errs  , r_hi_errs)

c = rt.TCanvas('c_fit', 'c_fit', 800, 800)
pad1 = rt.TPad('pad1', 'pad1', 0., 0.3, 1.0, 1.0); pad1.Draw()
pad2 = rt.TPad('pad2', 'pad2', 0., 0.0, 1.0, 0.3); pad2.Draw()
pad1.SetBottomMargin(0.01)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.25)

# Add the fit distribution
pad1.cd()
g_r.SetTitle("Z' fit rate vs. Z' mass;; #sigma(Z')#it{B}(Z'#rightarrowe#mu) (fb)")
g_r.SetMarkerStyle(20)
g_r.SetMarkerSize(0.8)
g_r.SetLineWidth(2)
g_r.SetMarkerColor(rt.kBlack)
g_r.SetLineColor(rt.kBlack)
g_r.Draw("APE1")

g_r.GetXaxis().SetRangeUser(xmin, xmax)
g_r.GetYaxis().SetRangeUser(min_r - 0.05*(max_r - min_r), max_r + 0.1*(max_r - min_r))
g_r.GetXaxis().SetLabelSize(0.)

# Add an approximate significance distribution below the fit results
pad2.cd()
significances = array('d')
sig_half = array('d')
sig_err  = array('d')
m_errs = array('d') #create buffers between mass points
for index in range(len(r_fits)):
   r_err = (r_lo_errs[index] if r_fits[index] > 0. else r_hi_errs[index])
   if r_err <= 0.:
      print "Fit errors for mass point %.1f are 0" % (masses[index])
      r_err = max([0.1, r_lo_errs[index], r_hi_errs[index]])
   significances.append((r_fits[index]) / r_err)
   sig_half.append(significances[-1]/2.)
   sig_err.append(significances[-1]/2.)
   if index < len(masses) - 1:
      m_errs.append(0.8*(masses[index+1] - masses[index])/2.)
   elif len(masses) > 1:
      m_errs.append(0.8*(masses[index] - masses[index-1])/2.)
   else:
      m_errs.append(1.)

g_sig = rt.TGraphErrors(len(masses), masses, sig_half, m_errs, sig_err)
g_sig.SetTitle(";Z' mass (GeV);#sigma_{r}")
g_sig.SetFillColor(rt.kAtlantic)
g_sig.Draw("AE2")
g_sig.GetXaxis().SetRangeUser(xmin, xmax)
g_sig.GetYaxis().SetRangeUser(min(-4, max(significances)), max(4, max(significances)))
g_sig.GetXaxis().SetLabelSize(0.08)
g_sig.GetYaxis().SetLabelSize(0.08)
g_sig.GetXaxis().SetTitleSize(0.1)
g_sig.GetXaxis().SetTitleOffset(0.9)
g_sig.GetYaxis().SetTitleSize(0.1)
g_sig.GetYaxis().SetTitleOffset(0.4)
line = rt.TLine(xmin, 0., xmax, 0.)
line.SetLineWidth(2)
line.SetLineStyle(rt.kDashed)
line.SetLineColor(rt.kBlack)
line.Draw('same')

c.SaveAs(figdir+'fits.png')


#----------------------------------------------
# Significance plot
#----------------------------------------------

g_sig   = rt.TGraph(len(masses), masses, r_sigs)

c = rt.TCanvas('c_sig', 'c_sig', 800, 600)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.08)
g_sig.SetTitle("Measurement significance vs. Z' mass; Z' mass (GeV); #sigma(#it{B}(Z'#rightarrowe#mu))")
g_sig.SetMarkerStyle(20)
g_sig.SetMarkerSize(0.75)
g_sig.SetLineWidth(2)
g_sig.SetMarkerColor(rt.kRed)
g_sig.SetLineColor  (rt.kGray+2)
g_sig.Draw("APL")

min_sig = min(r_sigs)
max_sig = max(r_sigs)

g_sig.GetXaxis().SetRangeUser(xmin, xmax)
g_sig.GetYaxis().SetRangeUser(min_sig - 0.1*(max_sig-min_sig), max_sig + 0.1*(max_sig-min_sig))

g_sig.GetXaxis().SetLabelSize(0.04)
g_sig.GetXaxis().SetTitleSize(0.045)
g_sig.GetXaxis().SetTitleOffset(1.00)
g_sig.GetYaxis().SetLabelSize(0.04)
g_sig.GetYaxis().SetTitleSize(0.045)
g_sig.GetYaxis().SetTitleOffset(0.65)

line = rt.TLine(xmin, 0., xmax, 0.)
line.SetLineColor(rt.kBlack)
line.SetLineStyle(rt.kDashed)
line.SetLineWidth(2)
line.Draw("same")

c.SaveAs(figdir+'sig.png')


#----------------------------------------------
# p-value plot
#----------------------------------------------

pvals = array('d')
for sig in r_sigs: pvals.append(rt.RooStats.SignificanceToPValue(sig))

# Calculate the global significance by evaluating N(up crossings)
ref_u = 0.
n_ref_u = 0

for index in range(1,len(r_sigs)):
   if r_sigs[index-1] < ref_u and r_sigs[index] > ref_u: n_ref_u += 1

if len(global_pvals) == 0:
   for pval in pvals:
      sig = rt.RooStats.PValueToSignificance(pval)
      global_pval = min(1., pval + n_ref_u*exp(-(sig**2-ref_u**2)/2.))
      global_pvals.append(global_pval)

g_pval   = rt.TGraph(len(masses), masses, pvals)
g_global = rt.TGraph(len(masses), masses, global_pvals)

c = rt.TCanvas('c_pval', 'c_pval', 800, 600)
g_pval.SetTitle("Measurement p-value vs. Z' mass; Z' mass (GeV); p")
g_pval.SetMarkerStyle(20)
g_pval.SetMarkerSize(0.8)
g_pval.SetLineWidth(2)
g_pval.SetMarkerColor(rt.kRed)
g_pval.SetLineColor(rt.kRed)
g_pval.Draw("AL")
# g_pval.Draw("APL")

g_global.SetMarkerStyle(20)
g_global.SetMarkerSize(0.8)
g_global.SetLineWidth(2)
g_global.SetMarkerColor(rt.kRed)
g_global.SetLineColor(rt.kRed)
g_global.SetLineStyle(rt.kDashed)
g_global.Draw("L")

min_pval = min(pvals)
max_pval = max(pvals)
g_pval.GetXaxis().SetRangeUser(xmin, xmax)
g_pval.GetYaxis().SetRangeUser(0.2*min_pval, 2.)
c.SetLogy()

# Add sigma lines to the p-value plot
sig_p_min = rt.RooStats.PValueToSignificance(0.2*min_pval)
lines = []
for sig in range(int(sig_p_min)+1):
   p_sig = rt.RooStats.SignificanceToPValue(sig)
   line = rt.TLine(xmin, p_sig, xmax, p_sig)
   line.SetLineStyle(rt.kDashed)
   line.SetLineWidth(2)
   line.SetLineColor(rt.kBlack)
   line.Draw("same")
   lines.append(line)

c.SaveAs(figdir+'pval.png')

# Make just a histogram of log(p-value)
h = rt.TH1D('h_pval', '-log_{10}(p-value) distribution', 40, 0, 8)
for pval in pvals: h.Fill(-log10(abs(pval)))
h.SetLineColor(rt.kBlue)
h.SetFillColor(rt.kAtlantic)
h.SetLineWidth(2)
h.Draw("hist")
h.SetXTitle('-log_{10}(p)')
c.SaveAs(figdir+'pval_hist.png')
