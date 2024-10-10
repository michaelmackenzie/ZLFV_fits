# Make PDF templates from an input histogram file
import os
import argparse
import ROOT as rt
from math import exp,pow,sqrt

#--------------------------------------------------------------------------------
# Double-sided Crystal Ball function
def dscb(x, mean, sigma, alpha1, enne1, alpha2, enne2):
  # function constants
  a1 = pow(abs(enne1/alpha1), enne1)*exp(-0.5*alpha1*alpha1)
  a2 = pow(abs(enne2/alpha2), enne2)*exp(-0.5*alpha2*alpha2)
  b1 = abs(enne1/alpha1) - abs(alpha1)
  b2 = abs(enne2/alpha2) - abs(alpha2)

  # function main variable
  var = (x - mean) /sigma
  val = 0.

  if  (var < -alpha1): val = a1*pow(b1 - var, -enne1)
  elif(var <      0.): val = exp(-0.5*var*var)
  elif(var <  alpha2): val = exp(-0.5*var*var) #same sigma for left/right
  else               : val = a2*pow(b2 + var, -enne2)

  return val

#--------------------------------------------------------------------------------
# Double-sided Crystal Ball call for TF1
def dscb_func(X, P):
  return P[0]*dscb(X[0], P[1], P[2], P[3], P[4], P[5], P[6]);

#--------------------------------------------------------------------------------
# Fit an input histogram
def fit_hist(hist, figdir, min_mass = 95., signal = 'zprime', set = 11):
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(1111)
    name = hist.GetName()
    if signal == 'zprime':
        if name == 'top':
            func = rt.TF1('func_'+name, 'pol2(0)', min_mass, 150)
        elif name == 'higgs' or name == 'zmumu':
            func = rt.TF1('func_'+name, 'TMath::Exp([0] + (x/100.)*[1]) + [2] + [3]*(x/100.) + [4]*TMath::Power(x/100.,2)', min_mass, 150)
            func.SetParameters(10., -10., 0.1, -1., 1.)
        else:
            func = rt.TF1('func_'+name, 'expo(0)', min_mass, 150)
    else:
        if name == 'top' or name == 'qcd' or name == 'higgs':
            func = rt.TF1('func_'+name, 'pol2(0)', 70., 110.)
        elif name == 'zmumu':
            func = rt.TF1("func_"+name, dscb_func, 70., 110., 7, 1);
            if set == 11:
                func.SetParameters(hist.Integral(), 80., 6., 1.4, 7.8, 1.9, 1.2);
            elif set == 12:
                func.SetParameters(hist.Integral(), 83., 5., 1.056, 7.490, 1.000, 7.597);
            elif set == 13:
                func.SetParameters(hist.Integral(), 83., 5., 1.056, 7.490, 1.000, 7.597);
            func.SetParNames("Norm", "#mu", "#sigma", "#alpha_{1}", "n_{1}", "#alpha_{2}", "n_{2}");
            func.SetParLimits(func.GetParNumber("#mu")       , 70., 90.);
            func.SetParLimits(func.GetParNumber("#sigma")    ,  3., 10.);
            func.SetParLimits(func.GetParNumber("n_{1}")     , 0.1, 10.);
            func.SetParLimits(func.GetParNumber("n_{2}")     , 0.1, 10.);
            func.SetParLimits(func.GetParNumber("#alpha_{1}"), 0.1,  5.);
            func.SetParLimits(func.GetParNumber("#alpha_{2}"), 0.1,  5.);
        else:
            func = rt.TF1("func_"+name, "[0]*TMath::Gaus(x, [1], [2], true) + [3] + [4]*x + [5]*x*x", 70., 110.);
            func.SetParameters(hist.Integral(), 65., 10., 70., -0.1, 0.);
            func.SetParLimits(0, 0., 1.e9);
            func.SetParLimits(1, 40., 70.);
            func.SetParLimits(2, 1., 30.);

    hist.Fit(func, 'R 0')
    c = rt.TCanvas()
    hist.SetLineWidth(2)
    hist.SetLineColor(rt.kBlue)
    hist.SetMarkerColor(rt.kBlue)
    hist.SetMarkerSize(0.8)
    hist.SetMarkerStyle(20)
    hist.Draw('E1')
    func.Draw('same')
    if name == 'ztautau':
        func_2 = func.Clone('func_2')
        func_2.SetParameter(0, 0.)
        func_2.SetLineStyle(rt.kDashed)
        func_2.Draw('same')
        
    hist.SetXTitle('M_{e#mu} (GeV/c^{2})')
    hist.SetYTitle('Events / %.1f GeV/c^{2}' % (hist.GetBinWidth(1)))
    hist.GetXaxis().SetRangeUser(min_mass, 150.)
    c.SaveAs(figdir+name+'_fit.png')
    hist.GetYaxis().SetRangeUser(1.e-3*hist.GetMaximum(), 5.*hist.GetMaximum())
    c.SetLogy()
    c.SaveAs(figdir+name+'_fit_log.png')
    return func

rt.gROOT.SetBatch(True)

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--set", dest="set",default=11, type=int,help="Histogram file set")
parser.add_argument("--tag", dest="tag",default="v01", type=str,help="Output naming tag")
parser.add_argument("--mass", dest="mass",default=110., type=float,help="Signal mass for the template")
parser.add_argument("--min-mass", dest="min_mass",default=95., type=float,help="Minimum mass used in the data card")
parser.add_argument("--max-mass", dest="max_mass",default=700., type=float,help="Maximum mass used in the data card")
parser.add_argument("--signal", dest="signal",default="zprime", type=str,help="Signal name")
parser.add_argument("--use-hists", dest="use_hists",default=False, action='store_true', help="Use the input histograms for the templates")
parser.add_argument("--use-pdfs", dest="use_pdfs",default=False, action='store_true', help="Use parametric PDFs in the output workspace")

args, unknown = parser.parse_known_args()
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()


#----------------------------------------------
# Output data paths
#----------------------------------------------

figdir = "./figures/templates/%s_set_%i_%s/" % (args.signal, args.set, args.tag)
outdir = "./templates/%s/" % (args.signal)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir))
os.system("[ ! -d %s ] && mkdir -p %s" % (outdir , outdir))

#----------------------------------------------
# Read in the input MC info
#----------------------------------------------

f_in = rt.TFile.Open('hists/%s_lepm_%i_2016_2017_2018.hist' % (args.signal, args.set), 'READ')
if args.signal == 'zemu':
    stack = f_in.Get('bkg_stack')
    top = stack.GetHists().At(0)
    ww = stack.GetHists().At(1)
    zmumu = stack.GetHists().At(2)
    ztautau = stack.GetHists().At(3)
    qcd = stack.GetHists().At(4)

    top.SetName('top')
    ww.SetName('ww')
    zmumu.SetName('zmumu')
    ztautau.SetName('ztautau')
    qcd.SetName('qcd')

    top.SetTitle('Top')
    ww.SetTitle('WW')
    zmumu.SetTitle('Z#rightarrow#mu#mu')
    ztautau.SetTitle('Z#rightarrow#tau#tau')
    qcd.SetTitle('QCD')

    hists = [ww, top, qcd, zmumu, ztautau]
else:
    stack = f_in.Get('hstack')
    higgs = stack.GetHists().At(0)
    top = stack.GetHists().At(1)
    ww = stack.GetHists().At(2)
    zmumu = stack.GetHists().At(3) # Not used here
    ztautau = stack.GetHists().At(4)
    qcd = stack.GetHists().At(5)

    higgs.SetName('higgs')
    top.SetName('top')
    ww.SetName('ww')
    zmumu.SetName('zmumu')
    ztautau.SetName('ztautau')
    qcd.SetName('qcd')

    higgs.SetTitle('H#rightarrow#tau#tau/WW')
    top.SetTitle('Top')
    ww.SetTitle('WW')
    zmumu.SetTitle('Z#rightarrow#mu#mu')
    ztautau.SetTitle('Z#rightarrow#tau#tau')
    qcd.SetTitle('QCD')

    if args.min_mass < 95.: #only include Z->mumu below this threshold
        hists = [ww, top, qcd, zmumu, ztautau, higgs]
    else:
        hists = [ww, top, qcd, ztautau, higgs]

#----------------------------------------------
# Process the histograms of interest
#----------------------------------------------

if not args.use_hists:
    fit_res = []
    for hist in hists: fit_res.append(fit_hist(hist, figdir, args.min_mass, args.signal, args.set))

#----------------------------------------------
# Make a background template from the fits
#----------------------------------------------

mass = args.mass
width = 1./50.*mass
if args.signal == 'zemu':
    min_mass = 70.
    max_mass = 110.
    nbins    = 80
else:
    min_mass = max(mass - 10.*width, args.min_mass)
    max_mass = min(mass + 10.*width, args.max_mass)
    nbins    = int((max_mass - min_mass) / (0.5 * width))
hbkg = rt.TH1D('hbkg', 'Background template', nbins, min_mass, max_mass)
h_no_zmm = rt.TH1D('h_no_zmm', 'Background template without Z->mumu', nbins, min_mass, max_mass)
h_zmm = rt.TH1D('h_zmm', 'Z->mumu template', nbins, min_mass, max_mass)
x_w_in = ww.GetBinWidth(1)
x_w_out = hbkg.GetBinWidth(1)
print x_w_in, x_w_out

for xbin in range(nbins):
    x = hbkg.GetBinCenter(xbin+1)
    bkg_val = 0.
    no_zmm_val = 0.
    if args.use_hists:
        for h in hists:
            bkg_val += max(0., h.GetBinContent(h.FindBin(x)))/x_w_in*x_w_out
            if h.GetName() != 'zmumu': no_zmm_val += max(0., h.GetBinContent(h.FindBin(x)))/x_w_in*x_w_out
    else:
        for fit in fit_res:
            # val = max(0., fit.Integral(x-x_w_in/2., x+x_w_in/2.)/x_w_in)*(x_w_out/x_w_in)
            val = max(0., fit.Eval(x))/x_w_in*x_w_out
            bkg_val += val
            if not 'zmumu' in fit.GetName(): no_zmm_val += val
    hbkg.SetBinContent(xbin+1, bkg_val)
    hbkg.SetBinError(xbin+1, sqrt(bkg_val))
    h_no_zmm.SetBinContent(xbin+1, no_zmm_val)
    h_no_zmm.SetBinError(xbin+1, sqrt(no_zmm_val))
    h_zmm.SetBinContent(xbin+1, bkg_val - no_zmm_val)
    h_zmm.SetBinError(xbin+1, sqrt(max(0., bkg_val - no_zmm_val)))

c = rt.TCanvas()
hbkg.Draw('hist')
hbkg.SetLineWidth(2)
hbkg.SetXTitle('M_{e#mu} (GeV/c^{2})')
hbkg.SetYTitle('Events / %.1f GeV/c^{2}' % (x_w_out))
c.SaveAs(figdir+'template.png')

#----------------------------------------------
# Make a data card for the MC template
#----------------------------------------------
cardname = 'datacard_' + args.signal + '_mc_template_'  + args.tag + '_' + str(args.set) + ('_mass-%.2f' % (mass)) + '.txt'
wsname   = 'workspace_' + args.signal + '_mc_template_' + args.tag + '_' + str(args.set) + ('_mass-%.2f' % (mass)) + '.root'

# Standard preamble
txt="# -*- mode: tcl -*-\n"
txt+="#Auto generated Z prime search COMBINE datacard\n"
txt+="#Using signal mass = %.3f\n\n" % (mass)
txt+="imax * number of channels\njmax * number of backgrounds\nkmax * number of nuisance parameters\n\n"

# Define the background and signal PDFs
txt+="#-----------------------------------------------------------------------------------------------------------\n"
txt+="shapes signal      bin1 %s ws:sig_pdf\n" % (wsname)
txt+="shapes background  bin1 %s ws:bkg_pdf\n" % (wsname)
txt+="shapes data_obs    bin1 %s ws:data_obs\n" % (wsname)
if args.use_pdfs:
    txt+="shapes zmumu       bin1 %s ws:zmm_pdf\n" % (wsname)
    

txt+="\n#-----------------------------------------------------------------------------------------------------------\n\n"
# Define the channel
txt+="bin         bin1\n"
txt+="observation  -1\n\n"
if args.use_pdfs:
    txt+="bin         bin1      bin1       bin1\n"
    txt+="process    signal   background  zmumu\n\n"
    txt+="process       0         1          2\n\n"
    txt+="rate          1         1          1\n" #Rate is taken from the _norm variables
else:
    txt+="bin         bin1      bin1\n"
    txt+="process    signal   background\n\n"
    txt+="process       0         1\n\n"
    txt+="rate          1         1\n" #Rate is taken from the _norm variables

# Define the uncertainties
txt+="#-----------------------------------------------------------------------------------------------------------\n"
if args.use_pdfs:
    txt+="SignalSys lnN 1.10      -          -\n"
else:
    txt+="SignalSys lnN 1.10      -\n"

# Write the file
with open(outdir+cardname,'w') as fl:
    fl.write(txt)
fl.close()

#----------------------------------------------
# Make a workspace for the datacard
#----------------------------------------------

ws = rt.RooWorkspace('ws')
obs = rt.RooRealVar('mass_ll' if args.signal == 'zprime' else 'lepm_%i' % (args.set), 'Dilepton mass', (min_mass+max_mass)/2., min_mass, max_mass)
obs.setBins(nbins)
bkg_data    = rt.RooDataHist('bkg_hist', 'Background template', obs, hbkg)
no_zmm_data = rt.RooDataHist('no_zmm_hist', 'Background template without Z->#mu#mu', obs, h_no_zmm)
zmm_data    = rt.RooDataHist('zmm_hist', 'Z->#mu#mu template', obs, h_zmm)

bkg_pdf     = rt.RooHistPdf('bkg_pdf'   , 'Background PDF'                  , obs, bkg_data   )

# Generate toy data (not used)
mc_data = bkg_pdf.generateBinned(obs, hbkg.Integral())
mc_data.SetName('data_obs')
ws.Import(mc_data)

if args.use_pdfs:
    # Fit the non-Z->mumu PDF
    b1    = rt.RooRealVar("b1", "b1", 0.06, -3., 3.)
    b2    = rt.RooRealVar("b2", "b2", 0.00, -1., 1.)
    b3    = rt.RooRealVar("b3", "b3", 0.00, -1., 1.)
    pol   = rt.RooPolynomial("pol", "pol",obs,rt.RooArgList(b1,b2));
    gs_mu = rt.RooRealVar("mu","mu",65.,50,69);
    gs_wd = rt.RooRealVar("wd","wd",8.4,5,20);
    gauss = rt.RooGaussian("gauss", "gauss",obs,gs_mu,gs_wd);
    ratio = rt.RooRealVar("ratio", "ratio",0.1,0,1);
    pdf   = rt.RooAddPdf("bkg_pdf", "bkg_pdf", rt.RooArgList(gauss,pol),rt.RooArgList(ratio))
    pdf.fitTo(no_zmm_data)

    # Fit the Z->mumu
    zmm_mu  = rt.RooRealVar ("zmm_mu", "zmm_mu", 80.)
    zmm_wd  = rt.RooRealVar ("zmm_wd", "zmm_wd",  6.)
    zmm_a1  = rt.RooRealVar ("zmm_a1", "zmm_a1", 1.4)
    zmm_a2  = rt.RooRealVar ("zmm_a2", "zmm_a2", 7.8)
    zmm_n1  = rt.RooRealVar ("zmm_n1", "zmm_n1", 1.9)
    zmm_n2  = rt.RooRealVar ("zmm_n2", "zmm_n2", 1.2)
    zmm_pdf = rt.RooDoubleCB("zmm_pdf", "zmm_pdf", obs, zmm_mu, zmm_wd, zmm_a1, zmm_n1, zmm_a2, zmm_n2)
    zmm_pdf.fitTo(zmm_data)

    no_zmm_norm = rt.RooRealVar('bkg_pdf_norm', 'Background normalization', h_no_zmm.Integral(), 0., 1.e6)
    zmm_norm = rt.RooRealVar('zmm_pdf_norm', 'Z->mumu normalization', h_zmm.Integral(), 0., 1.e4)
    ws.Import(pdf)
    ws.Import(no_zmm_norm)
    ws.Import(zmm_pdf)
    ws.Import(zmm_norm)
else:
    bkg_norm = rt.RooRealVar('bkg_pdf_norm', 'Background normalization', hbkg.Integral(), 0., 1.e6)
    ws.Import(bkg_pdf)
    ws.Import(bkg_norm)

# Create an example signal PDF
sig_mean  = rt.RooRealVar('sig_mean', 'Signal mean', mass)
sig_width = rt.RooRealVar('sig_width', 'Signal width', mass/50.)
sig_mean.setConstant(True)
sig_width.setConstant(True)
sig_pdf = rt.RooGaussian('sig_pdf', 'Signal PDF', obs, sig_mean, sig_width)
sig_norm = rt.RooRealVar('sig_pdf_norm', 'Background normalization', 50.)
sig_norm.setConstant(True)
ws.Import(sig_pdf)
ws.Import(sig_norm)
ws.writeToFile(outdir+wsname)


# Test fitting the non-Z->mumu template with a Gaussian + polynomial

b1    = rt.RooRealVar("b1", "b1", 0.06, -3., 3.)
b2    = rt.RooRealVar("b2", "b2", 0.00, -1., 1.)
b3    = rt.RooRealVar("b3", "b3", 0.00, -1., 1.)
pol   = rt.RooPolynomial("pol", "pol",obs,rt.RooArgList(b1,b2));
gs_mu = rt.RooRealVar("mu","mu",65.,50,69);
gs_wd = rt.RooRealVar("wd","wd",8.4,5,20);
gauss = rt.RooGaussian("gauss", "gauss",obs,gs_mu,gs_wd);
ratio = rt.RooRealVar("ratio", "ratio",0.1,0,1);
pdf   = rt.RooAddPdf("pdf", "pdf", rt.RooArgList(gauss,pol),rt.RooArgList(ratio))

# Add a signal PDF
sig_mu = rt.RooRealVar("sig_mu", "sig_mu", 90.69); sig_mu.setConstant(True)
sig_wd = rt.RooRealVar("sig_wd", "sig_wd",  2.32); sig_wd.setConstant(True)
sig_a1 = rt.RooRealVar("sig_a1", "sig_a1",  1.01); sig_a1.setConstant(True)
sig_a2 = rt.RooRealVar("sig_a2", "sig_a2",  1.34); sig_a2.setConstant(True)
sig_n1 = rt.RooRealVar("sig_n1", "sig_n1",  3.08); sig_n1.setConstant(True)
sig_n2 = rt.RooRealVar("sig_n2", "sig_n2",  2.57); sig_n2.setConstant(True)
sig_pdf = rt.RooDoubleCB("sig_pdf", "sig_pdf", obs, sig_mu, sig_wd, sig_a1, sig_n1, sig_a2, sig_n2)

# Combine the signal and background
nbkg = rt.RooRealVar("nbkg", "nbkg", 1.e4, 0., 1.e6)
nsig = rt.RooRealVar("nsig", "nsig", 0., -1.e3, 1.e3) #; nsig.setConstant(True)
tot_pdf = rt.RooAddPdf("tot_pdf", "Total PDF", rt.RooArgList(pdf, sig_pdf), rt.RooArgList(nbkg, nsig))

# Fit the non-Z->mumu template
no_zmm_pdf  = rt.RooHistPdf('no_zmm_pdf', 'Background PDF without Z->#mu#mu', obs, no_zmm_data)
zmm_pdf     = rt.RooHistPdf('zmm_pdf'   , 'Z->#mu#mu PDF'                   , obs, zmm_data   )
rt.Math.MinimizerOptions.SetDefaultPrecision(1e-15)
tot_pdf.fitTo(no_zmm_data)

# Plot the results
frame = obs.frame()
no_zmm_data.plotOn(frame, rt.RooFit.Name("data"))
tot_pdf.plotOn(frame, rt.RooFit.LineColor(rt.kRed), rt.RooFit.Name("pdf"))
tot_pdf.plotOn(frame, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.LineStyle(rt.kDashed), rt.RooFit.Components("pol"))
tot_pdf.plotOn(frame, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.LineStyle(rt.kDashed), rt.RooFit.Components("gauss"))
# no_zmm_pdf.plotOn(frame, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.Name("template"))

c = rt.TCanvas("c","c", 700, 800)
pad1 = rt.TPad("pad1", "pad1", 0., 0.3, 1.0, 1.0); pad1.Draw()
pad2 = rt.TPad("pad2", "pad2", 0., 0.0, 1.0, 0.3); pad2.Draw()

pad1.cd()
frame.Draw()
pad2.cd()
hpull = frame.pullHist("data","pdf");
hpull.SetTitle("")
hpull.Draw()
hpull.GetXaxis().SetRangeUser(args.min_mass, args.max_mass-1.e-3)
hpull.GetYaxis().SetRangeUser(-2,2.)
hpull.GetXaxis().SetLabelSize(0.08)
hpull.GetYaxis().SetLabelSize(0.08)
hpull.SetLineColor(rt.kRed)
# hpull.SetMarkerColor(rt.kRed)

# hpull_2 = frame.pullHist("data","template");
# hpull_2.SetLineColor(rt.kBlue)
# hpull_2.SetMarkerColor(rt.kBlue)
# hpull_2.SetMarkerStyle(4)
# hpull_2.Draw("PE1 same")

c.SaveAs(figdir+'template_fit.png')
pad1.SetLogy()
c.SaveAs(figdir+'template_fit_log.png')


print "N(bkg) = %.1f, N(sig) = %.3f" % (nbkg.getVal(), nsig.getVal())
