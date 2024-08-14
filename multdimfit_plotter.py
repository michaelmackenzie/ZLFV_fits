import ROOT
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-f', dest='flname', default="",type=str, help='sum the integers (default: find the max)')
parser.add_argument('-o', dest='output', default="",type=str, help='sum the integers (default: find the max)')
parser.add_argument('--mllvar', dest='mllvar', default="",type=str, help='sum the integers (default: find the max)')
parser.add_argument('--bin', dest='bin_name', default="",type=str, help='sum the integers (default: find the max)')
parser.add_argument('-v', dest='verb', default=False, action='store_true', help='sum the integers (default: find the max)')
args = parser.parse_args()



f = ROOT.TFile('/afs/cern.ch/work/g/gkaratha/private/Analysis/DispJets/Analyzer/Limit/TutorialCombine/CMSSW_11_3_4/src/combinetutorial-2023-parametric/'+args.flname)
w = f.Get("w")
if args.verb:
   w.Print("v")

n_bins = 80
binning = ROOT.RooFit.Binning(n_bins,70,110)

can = ROOT.TCanvas()
plot = w.var(args.mllvar).frame()
w.data("data_obs").plotOn( plot, binning )

# Load the S+B model
sb_model = w.pdf("model_s").getPdf(args.bin_name)

# Prefit
sb_model.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit") )

# Postfit
w.loadSnapshot("MultiDimFit")
sb_model.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit") )
r_bestfit = int(w.var("r").getVal()*100)/100.0
r_error = int(w.var("r").getError()*100)/100.0
print ("r=",r_bestfit,"+/-",r_error)
plot.Draw()

leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
leg.AddEntry("prefit", "Prefit S+B model (r=1.00)", "L")
leg.AddEntry("postfit", "Postfit S+B model (r={0} +/- {1})".format(r_bestfit,r_error), "L")
leg.Draw("Same")

can.Update()
can.SaveAs("combine_fit_"+args.output+".png")
