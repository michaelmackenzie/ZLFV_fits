# Plot the Z prime BDT score for various signal masses
import os
import argparse
import ROOT as rt
from array import array

def cdf_from_bdt(bdt, name):
    cdf = bdt.Clone('cdf_'+name)
    nbins = bdt.GetNbinsX()
    for index in range(nbins):
        integral = bdt.Integral(index+1, nbins)
        cdf.SetBinContent(index+1, integral)
        cdf.SetBinError(index+1, 0.)
    return cdf

rt.gROOT.SetBatch(True)

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/signal_bdt/"
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

### MC signal mass points
sgn_masspoints=["100", "200","400","600","800","1000"]
sgn_masspoint_files={
               "100":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun18.root",\
               "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
               "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
               "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}
ndens={"100":99200.,"200":96300.,"400":97600.,"600":97600.,"800":97600.,"1000":97600.}

# Loop through the signal files
bdts = []
cdfs = []
for mpoint in sgn_masspoints:
    # Create a TChain for the signal
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_masspoint_files[mpoint])

    # Create a high binning of the BDT score for CDF creation
    hname = "hbdt_"+mpoint
    bdt = rt.TH1F(hname,"BDT score",3000,0,1)
    cc.Draw("xgb>>"+hname,"(Flag_met && Flag_muon && Flag_electron)")
    
    scale = 1. / ndens[mpoint] # divide out N(gen)
    bdt.Scale(scale)
    efficiency = bdt.Integral() # Acceptance before BDT score cuts
    bdt.Scale(1./efficiency) # Remove acceptance efficiency
    print "Mass point %s: Efficiency = %.3f" % (mpoint, efficiency)
    bdt.SetTitle("Mass = %s GeV" % (mpoint))

    cdfs.append(cdf_from_bdt(bdt, mpoint))
    bdt.Rebin(100)
    bdts.append(bdt)

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

c.SaveAs(figdir+'bdt.png')

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

c.SaveAs(figdir+'cdf.png')
