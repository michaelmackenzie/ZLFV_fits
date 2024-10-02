# Optimize the BDT score category as a function of mass
import os
import argparse
import ROOT as rt
from array import array
from math import sqrt

#----------------------------------------------------------------------------
# Interpolate the CDF linearly in bin height
def interpolate_cdf(mass, cdf_1, mass_1, cdf_2, mass_2):
    cdf   = cdf_1.Clone('cdf')
    nbins = int(cdf.GetNbinsX())
    for index in range(1,nbins+1):
        val_1 = cdf_1.GetBinContent(index)
        val_2 = cdf_2.GetBinContent(index)
        if abs(mass_1 - mass_2) < 1.: val = val_1
        else: val = val_1 + (val_2 - val_1)/(mass_2 - mass_1)*(mass - mass_1)
        val = max(0., min(1., val))
        cdf.SetBinContent(index, val)
    return cdf

#----------------------------------------------------------------------------
# BDT score --> CDF transform
def cdf_from_bdt(bdt, name):
    cdf   = bdt.Clone('cdf_'+name)
    nbins = int(bdt.GetNbinsX())
    for index in range(nbins):
        integral = bdt.Integral(index+1, nbins)
        cdf.SetBinContent(index+1, integral)
        cdf.SetBinError(index+1, 0.)
    return cdf


# Optimization outline:
# - Assume the background is dominantly WW/ttbar->WW
# - Need N(background events) and eff(signal events) as a function of mass and BDT score
# - Define the BDT score cut to optimize S/sqrt(B+S) ~ S/sqrt(B)
# - Make a distribution of N(data events) vs mass and CDF(data) vs. mass, use independently

#----------------------------------------------
# Initialize paths and datasets
#----------------------------------------------

rt.gROOT.SetBatch(True)

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/bdt_bin_optimization/"
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir))

### MC signal mass points
sgn_masspoints=["200","400","600","800","1000"]
sgn_masspoint_files={
               "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
               "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
               "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}
ndens={"200":96300.,"400":97600.,"600":97600.,"800":97600.,"1000":97600.}

### Data
data_file="Meas_full_bdt_v7_emu_scan_data_Run1*.root"
lumi=137.6

#----------------------------------------------
# Read in the signal info
#----------------------------------------------

# Loop through the signal files
bdts = []
cdfs = []
cuts = "(Flag_met && Flag_muon && Flag_electron)"
for mpoint in sgn_masspoints:
    # Create a TChain for the signal
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_masspoint_files[mpoint])

    # Create a high binning of the BDT score for CDF creation
    hname = "hbdt_"+mpoint
    bdt = rt.TH1F(hname,"BDT score",3000,0,1)
    cc.Draw("xgb>>"+hname,cuts)
    
    scale = 1. / ndens[mpoint] # divide out N(gen)
    bdt.Scale(scale)
    efficiency = bdt.Integral() # Acceptance before BDT score cuts
    bdt.Scale(1./efficiency) # Remove acceptance efficiency
    print "Mass point %s: Efficiency = %.3f" % (mpoint, efficiency)
    bdt.SetTitle("Mass = %s GeV" % (mpoint))

    cdfs.append(cdf_from_bdt(bdt, mpoint))
    bdt.Rebin(100)
    bdts.append(bdt)

#----------------------------------------------
# Read in the data
#----------------------------------------------

# Loop through the relevant mass points for data
data_cdfs = []
# Create a TChain for the data
cc=rt.TChain("mytreefit")
cc.Add(path+"/"+data_file)
data_mass = rt.TH1F("hdata_mass","Mass distribution",100, 100, 1000)
cc.Draw("mass_ll>>hdata_mass",cuts)
data_masspoints = ["200", "300", "400", "500"]
for mpoint in data_masspoints:
    mass = float(mpoint)
    mass_range = 50.

    # Create a high binning of the BDT score for CDF creation
    hname = "hdata_bdt_"+mpoint
    bdt = rt.TH1F(hname,"BDT score",3000,0,1)
    cc.Draw("xgb>>"+hname,cuts + " && mass_ll > " + str(mass-mass_range) + " && mass_ll < " + str(mass+mass_range))
    
    bdt.Scale(1./bdt.Integral())
    print "Mass point %s: Efficiency = %.3f" % (mpoint, efficiency)
    bdt.SetTitle("Mass = %s GeV" % (mpoint))
    cdf = cdf_from_bdt(bdt, "data_"+mpoint)
    data_cdfs.append(cdf)

#--------------------------------------------------------------------
# Plot the BDT scores
#--------------------------------------------------------------------

c = rt.TCanvas('c_bdt', 'c_bdt', 700, 500)
haxis = bdts[0]
colors = [rt.kBlue, rt.kGreen, rt.kRed, rt.kOrange, rt.kBlack, rt.kMagenta, rt.kAtlantic]
leg = rt.TLegend(0.1, 0.75, 0.9, 0.9)
leg.SetTextSize(0.04)
leg.SetNColumns(3)
for index in range(len(bdts)):
    h = bdts[index]
    if index == 0: h.Draw("hist")
    else         : h.Draw("hist same")
    h.SetLineColor(colors[index % len(colors)])
    h.SetLineWidth(2)
    leg.AddEntry(h)
leg.Draw()

haxis.SetTitle("BDT score vs. Mass")
haxis.SetXTitle("BDT score")
haxis.GetYaxis().SetRangeUser(0., 1.3*haxis.GetMaximum())
rt.gStyle.SetOptStat(0)

c.SaveAs(figdir+'signal_bdt.png')

#--------------------------------------------------------------------
# Plot the BDT score CDFs
#--------------------------------------------------------------------

c = rt.TCanvas('c_cdf', 'c_cdf', 700, 500)
haxis = cdfs[0]
for index in range(len(cdfs)):
    h = cdfs[index]
    if index == 0: h.Draw("hist")
    else         : h.Draw("hist same")
    h.SetLineColor(colors[index % len(colors)])
    h.SetLineWidth(2)
leg.Draw()

haxis.SetTitle("BDT score CDF vs. Mass")
haxis.SetXTitle("BDT score")
haxis.SetYTitle("CDF")
haxis.GetYaxis().SetRangeUser(0., 1.3*haxis.GetMaximum())
rt.gStyle.SetOptStat(0)

c.SaveAs(figdir+'signal_cdf.png')

#--------------------------------------------------------------------
# Plot the data CDFs
#--------------------------------------------------------------------

leg = rt.TLegend(0.1, 0.75, 0.9, 0.9)
leg.SetTextSize(0.04)
leg.SetNColumns(3)
haxis = data_cdfs[0]
for index in range(len(data_cdfs)):
    h = data_cdfs[index]
    if index == 0: h.Draw("hist")
    else         : h.Draw("hist same")
    h.SetLineColor(colors[index % len(colors)])
    h.SetLineWidth(2)
    leg.AddEntry(h)
leg.Draw()

haxis.SetTitle("Data BDT score CDF vs. Mass")
haxis.SetXTitle("BDT score")
haxis.SetYTitle("CDF")
haxis.GetYaxis().SetRangeUser(0., 1.3*haxis.GetMaximum())
rt.gStyle.SetOptStat(0)

c.SaveAs(figdir+'data_cdf.png')

#--------------------------------------------------------------------
# Plot the data mass
#--------------------------------------------------------------------

func = rt.TF1("func", "exp([0] + [1]*x) + exp([2] + [3]*x) + exp([4] + [5]*x)", 100., 1000.)
# func.SetParameters(1., -0.01, 1, -0.04, 1., -0.001)
func.SetParameters(9., -0.0167, 11., -0.0167, 8., -0.008)
data_mass.Fit(func, "R")

data_mass.Draw("E1")
data_mass.SetMarkerStyle(20)
data_mass.SetMarkerSize(0.8)
data_mass.SetLineWidth(2)
data_mass.SetTitle("Data Mass distribution (no BDT cut)")
data_mass.SetXTitle("Mass (GeV/c^{2})")
data_mass.SetYTitle("")
data_mass.GetYaxis().SetRangeUser(0., 1.1*data_mass.GetMaximum())
func.Draw("same")
rt.gStyle.SetOptStat(0)

c.SaveAs(figdir+'data_mass.png')


#--------------------------------------------------------------------
# Test a single optimization point
#--------------------------------------------------------------------

mass = 200.
cdf_sig = cdfs[0]
cdf_data = data_cdfs[0]
rate_data = func.Eval(mass) / data_mass.GetBinWidth(1) * (2.*4./200.*mass) # Background rate in a +-1 sigma window around the signal peak
max_steps = 20
bdt_cuts = array('d')
sigs = array('d')
max_val = -1.
best_cut = 0.
for step in range(max_steps):
    bdt_cut = step * 1./max_steps
    sig_eff = cdf_sig.GetBinContent(cdf_sig.FindBin(bdt_cut))
    data_eff = cdf_data.GetBinContent(cdf_data.FindBin(bdt_cut))*rate_data
    ndata = max(3., rate_data * data_eff)
    sig = sig_eff / sqrt(ndata)
    sigs.append(sig)
    bdt_cuts.append(bdt_cut)
    if max_val < sig:
        max_val = sig
        best_cut = bdt_cut

g = rt.TGraph(len(bdt_cuts), bdt_cuts, sigs)
g.Draw("APL")
g.SetTitle("BDT cut optimization for %.0f GeV mass;BDT score cut;S/#sqrt{B}" % mass)
g.SetLineWidth(2)
g.SetLineColor(rt.kBlue)
g.SetMarkerSize(0.8)
g.SetMarkerStyle(20)
g.GetYaxis().SetRangeUser(0., 1.2*max_val)

c.SaveAs(figdir+'test_optimization.png')

#--------------------------------------------------------------------
# Test many optimization points
#--------------------------------------------------------------------

best_cuts = array('d')
masses = array('d')
min_mass = 200.
max_mass = 500.
max_steps = 50
subdir = figdir + 'optimizations/'
os.system("[ ! -d %s ] && mkdir -p %s" % (subdir , subdir))
for step in range(max_steps+1):
    mass = min_mass + (max_mass - min_mass)/max_steps*step

    # Find nearest signal CDF
    distance = 1.e10
    closest = sgn_masspoints[0]
    cl_index = 0
    for index in range(len(sgn_masspoints)):
        mpoint = sgn_masspoints[index]
        if abs(int(mpoint) - mass) < distance:
            distance = abs(int(mpoint) - mass)
            closest = mpoint
            cl_index = index
    cdf_1 = cdfs[cl_index]
    cdf_2_index = min(len(cdfs)-1, max(0, cl_index+1 if int(closest) < mass else cl_index-1))
    cdf_sig = interpolate_cdf(mass, cdf_1, float(closest), cdfs[cdf_2_index], float(sgn_masspoints[cdf_2_index]))
    cdf_sig.SetName("signal_cdf")

    # Find nearest data CDF
    distance = 1.e10
    closest = data_masspoints[0]
    cl_index = 0
    for index in range(len(data_masspoints)):
        mpoint = data_masspoints[index]
        if abs(int(mpoint) - mass) < distance:
            distance = abs(int(mpoint) - mass)
            closest = mpoint
            cl_index = index
    cdf_1 = data_cdfs[cl_index]
    cdf_2_index = min(len(data_cdfs)-1, max(0, cl_index+1 if int(closest) < mass else cl_index-1))
    cdf_data = interpolate_cdf(mass, cdf_1, float(closest), data_cdfs[cdf_2_index], float(data_masspoints[cdf_2_index]))
    cdf_data.SetName("data_cdf")

    cdf_sig.SetLineColor(rt.kBlue)
    cdf_sig.SetLineWidth(2)
    cdf_data.SetLineColor(rt.kRed)
    cdf_data.SetLineWidth(2)
    cdf_sig.Draw('L')
    cdf_data.Draw('L same')
    cdf_sig.GetYaxis().SetRangeUser(0.,1.1)
    c.SaveAs(subdir+'cdfs_%i.png' % (step))
    

    # Get the data rate for the optimization
    rate_data = func.Eval(mass) / data_mass.GetBinWidth(1) * (2.*4./200.*mass) # Background rate in a +-1 sigma window around the signal peak
    max_bdt_steps = 100
    max_val = -1.
    best_cut = 0.
    bdt_cuts = array('d')
    sigs = array('d')
    for bdt_step in range(max_bdt_steps):
        bdt_cut = bdt_step * 1./max_bdt_steps
        sig_eff = cdf_sig.GetBinContent(cdf_sig.FindBin(bdt_cut))
        data_eff = cdf_data.GetBinContent(cdf_data.FindBin(bdt_cut))*rate_data
        ndata = max(3., rate_data * data_eff)
        sig = sig_eff / sqrt(ndata)
        sigs.append(sig)
        bdt_cuts.append(bdt_cut)
        if max_val < sig:
            max_val = sig
            best_cut = bdt_cut

    best_cuts.append(best_cut)
    masses.append(mass)

    g = rt.TGraph(len(bdt_cuts), bdt_cuts, sigs)
    g.Draw("APL")
    g.SetTitle("BDT cut optimization for %.0f GeV mass;BDT score cut;S/#sqrt{B}" % mass)
    g.SetLineWidth(2)
    g.SetLineColor(rt.kBlue)
    g.SetMarkerSize(0.8)
    g.SetMarkerStyle(20)
    g.GetYaxis().SetRangeUser(0., 1.1*max_val)
    c.SaveAs(subdir+'optimization_%i.png' % (step))

g = rt.TGraph(len(masses), masses, best_cuts)
g.Draw("APL")
g.SetTitle("BDT cut optimization by mass;Mass (GeV/c^{2});BDT score cut")
g.SetLineWidth(2)
g.SetLineColor(rt.kBlue)
g.SetMarkerSize(0.8)
g.SetMarkerStyle(20)
g.GetYaxis().SetRangeUser(0., 1.)

c.SaveAs(figdir+'cut_optimization.png')

