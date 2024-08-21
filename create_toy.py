# Create a toy data set for Z prime emu resonance search testing
import os
import argparse
import ROOT as rt
from array import array

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o"           , dest="name"       , default="toy" , type=str  , help="Name for output")
parser.add_argument("--fit-file"   , dest="fit_file"   , default=""    , type=str  , help="File with input fit for toy generation")
parser.add_argument("--param"      , dest="param"      , default="bin1", type=str  , help="Parameter base name")
parser.add_argument("--seed"       , dest="seed"       , default=90    , type=int  , help="Toy generation seed")
parser.add_argument("--toy"        , dest="toy"        , default=1     , type=int  , help="Toy number")
parser.add_argument("--index"      , dest="index"      , default=-1    , type=int  , help="RooMultPdf index, -1 to use the default")
parser.add_argument("--signal-rate", dest="signal_rate", default=0.    , type=float, help="Signal rate to inject")
parser.add_argument("--signal-mass", dest="signal_mass", default=150.  , type=float, help="Signal mass to inject")

args, unknown = parser.parse_known_args()

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()

figdir = "./figures/toys_%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
   

rt.gROOT.SetBatch(True)

#----------------------------------------------
# Get the input background model
#----------------------------------------------

# Open the fit file
f = rt.TFile.Open(args.fit_file, 'READ')
ws = f.Get('ws_bkg')
pdf = ws.pdf('multipdf_%s' % (args.param))
if args.index >= 0:
   cat = ws.obj('pdfindex_%s' % (args.param))
   cat.setIndex(args.index)
mass = ws.var('mass_ll')
norm = ws.var('multipdf_%s_norm' % (args.param))

n_bkg = norm.getVal()

print "Nominal background rate = %.1f" % (n_bkg)


#----------------------------------------------
# Generate a toy dataset
#----------------------------------------------

rt.RooRandom.randomGenerator().SetSeed(args.seed)
dataset = pdf.generate(mass, rt.RooFit.NumEvents(rt.RooRandom.randomGenerator().Poisson(n_bkg)))
dataset.SetName("data_obs")

# Inject a signal if requested
if args.signal_rate > 0.:
   mean  = rt.RooRealVar("mean" , "Gaussian mean" , args.signal_mass)
   width = rt.RooRealVar("width", "Gaussian width", args.signal_mass/200.*4.)
   mean.setConstant(True)
   width.setConstant(True)
   sig_pdf = rt.RooGaussian("sig_pdf", "Signal PDF", mass, mean, width)
   if args.param == 'bin1':
      eff = max(0.1, 0.25 + 0.04*(args.signal_mass - 200.)/(800.))
   else:
      eff = max(0.1, 0.44 - 0.08*(args.signal_mass - 200.)/(800.))
   n_ref = 138. #N(signal) before efficiency cuts
   n_sig = args.signal_rate * n_ref * eff
   sig_data = sig_pdf.generate(mass, int(n_sig))
   dataset.append(sig_data)

frame = mass.frame(rt.RooFit.Title("Toy dataset %i" % (args.toy)))
dataset.plotOn(frame)
c = rt.TCanvas()
frame.Draw()
c.SaveAs(figdir+'toy_data_%i.png' % (args.toy))

os.system('[ ! -d WorkspaceScanTOY ] && mkdir -p WorkspaceScanTOY')
ws_out = rt.RooWorkspace('ws_bkg', 'Toy workspace')
ws_out.Import(dataset)
ws_out.writeToFile('WorkspaceScanTOY/toy_file_%s_%i.root' % (args.name, args.toy))
