# Plot the Z prime signal model interpolation
import os
import argparse
import ROOT as rt
from array import array

rt.gROOT.SetBatch(True)

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/signal_model/"
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

### MC signal mass points
sgn_masspoints=["90","200","400","600","800","1000"]
sgn_masspoint_files={
               "90" :"../NewElectronID/Meas_fullAndSF_bdt_v7_signal_mcRun18_updateID.root",\
               "200":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "400":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "600":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root",\
               "800":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root",\
               "1000":"Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root"}
# ndens={"90":208757,"200":44844.,"400":55294.,"600":60192.,"800":64756.,"1000":67263.}
#approximate values assuming gen sizes
ndens={"90":202933,"200":90856.,"400":110740.,"600":113408,"800":200000.,"1000":187676}

# Loop through the signal files
eff_array    = array('d')
mpoint_array = array('d')
hmasses = []
for mpoint in sgn_masspoints:
    # Create a TChain for the signal
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_masspoint_files[mpoint])

    # Create a high binning of the BDT score for CDF creation
    hname = "hmass_"+mpoint
    h = rt.TH1F(hname,"Di-lepton mass",1500,0,1500)
    cc.Draw("mass_ll>>"+hname,"(Flag_met && Flag_muon && Flag_electron)")
    
    scale = 1. / ndens[mpoint] # divide out N(gen)
    h.Scale(scale)
    h.SetTitle("Mass = %s GeV" % (mpoint))

    mpoint_array.append(float(mpoint))
    eff_array.append(h.Integral())
    hmasses.append(h)

print mpoint_array
print eff_array
#--------------------------------------------------------------------
# Plot the efficiencies
#--------------------------------------------------------------------

c = rt.TCanvas('c_eff', 'c_eff', 700, 500)
g = rt.TGraph(len(sgn_masspoints), mpoint_array, eff_array)
g.SetTitle("Total efficiency vs. mass; Z prime mass (GeV/c^{2});Efficiency")
g.SetMarkerStyle(20)
g.SetMarkerSize(0.8)
g.Draw("AP")

func = rt.TF1("func", "[0] + [1]*(x/1000) + [2]*pow(x/1000.,2)", 0., 1000.)
func.SetParameters(0.5, 0., 0.)
g.Fit(func, "R")
func.Draw("same")
c.SaveAs(figdir+"eff.png")




#--------------------------------------------------------------------
# Plot the mass distributions along with fits to them
#--------------------------------------------------------------------

fit_results = []
fit_errs = []
mass_errs = array('d')
for h in hmasses:
    mass_errs.append(0.)
    # Create a RooFit dataset and corresponding PDF to fit the signal distribution
    h_mean = h.GetMean()
    h_width = h.GetStdDev()
    min_mass = h_mean - 5.*h_width if h_mean > 100. else 70.
    max_mass = h_mean + 5.*h_width if h_mean > 100. else 110.
    obs = rt.RooRealVar("obs", "obs", h_mean, min_mass, max_mass, "GeV/c^{2}")
    dh = rt.RooDataHist("dh", "Mass distribution", obs, h)

    mean   = rt.RooRealVar("mean", "mean", h_mean, 0.9*h_mean, 1.1*h_mean)
    sigma  = rt.RooRealVar("sigma", "sigma", h_mean/50., h_mean/100., h_mean/25.);
    alpha1 = rt.RooRealVar("alpha1", "alpha1", 1., 0.1, 10.);
    alpha2 = rt.RooRealVar("alpha2", "alpha2", 1., 0.1, 10.);
    enne1  = rt.RooRealVar("enne1", "enne1", 5., 0.1, 30.);
    enne2  = rt.RooRealVar("enne2", "enne2", 5., 0.1, 30.);
    pdf    = rt.RooDoubleCrystalBall("pdf", "PDF", obs, mean, sigma, alpha1, enne1, alpha2, enne2)
    pdf.fitTo(dh, rt.RooFit.PrintLevel(-1), rt.RooFit.Warnings(0), rt.RooFit.PrintEvalErrors(-1), rt.RooFit.SumW2Error(1))

    # Print the fit results
    frame = obs.frame(rt.RooFit.Title("Signal model fit"))
    dh.plotOn(frame)
    pdf.plotOn(frame)

    c = rt.TCanvas()
    frame.Draw()
    c.SaveAs(figdir+h.GetName()+"_fit.png")

    fit_results.append([mean.getVal(), sigma.getVal(), alpha1.getVal(), alpha2.getVal(), enne1.getVal(), enne2.getVal()])
    fit_errs   .append([mean.getError(), sigma.getError(), alpha1.getError(), alpha2.getError(), enne1.getError(), enne2.getError()])

param_names = ["#mu", "#sigma", "#alpha_{1}", "#alpha_{2}", "n_{1}", "n_{2}"]

for param in range(len(param_names)):
    name = param_names[param]
    vals = array('d')
    errs = array('d')
    for index in range(len(fit_results)):
        vals.append(fit_results[index][param])
        errs.append(fit_errs   [index][param])
    c = rt.TCanvas('c_param', 'c_param', 700, 500)
    g = rt.TGraphErrors(len(vals), mpoint_array, vals, mass_errs, errs)
    g.SetTitle("Model %s parameter vs. mass;Z prime mass (GeV/c^{2});%s" % (name, name))
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.8)
    g.SetLineWidth(2)
    g.Draw("APE")
    min_val = min([vals[index]-errs[index] for index in range(len(vals))])
    max_val = max([vals[index]+errs[index] for index in range(len(vals))])
    g.GetYaxis().SetRangeUser(min_val - 0.1*(max_val-min_val), max_val + 0.1*(max_val-min_val))

    func = rt.TF1("func", "[0] + [1]*(x/1000) + [2]*pow(x/1000.,2)", 0., 1000.)
    func.SetParameters(0.5, 0., 0.)
    g.Fit(func, "R")
    func.Draw("same")
    c.SaveAs("%sparam_%i.png" % (figdir, param))



# haxis = bdts[0]
# colors = [rt.kBlue, rt.kGreen, rt.kRed, rt.kOrange, rt.kBlack, rt.kMagenta, rt.kAtlantic]
# leg = rt.TLegend(0.1, 0.75, 0.9, 0.9)
# leg.SetTextSize(0.04)
# leg.SetNColumns(3)
# for index in range(len(bdts)):
#     h = bdts[index]
#     if index == 0: h.Draw("hist")
#     else         : h.Draw("hist same")
#     h.SetLineColor(colors[index % len(colors)])
#     h.SetLineWidth(2)
#     leg.AddEntry(h)
# leg.Draw()

# haxis.SetTitle("BDT score vs. Mass")
# haxis.SetXTitle("BDT score")
# haxis.GetYaxis().SetRangeUser(0., 1.3*haxis.GetMaximum())
# rt.gStyle.SetOptStat(0)

# c.SaveAs(figdir+'bdt.png')
