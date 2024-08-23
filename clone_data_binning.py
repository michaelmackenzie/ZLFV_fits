# Clone the data binning for a toy data set
import os
import argparse
import ROOT as rt
from array import array

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--data", dest="data",default="", type=str,help="Input data file")
parser.add_argument("--toy", dest="toy",default="", type=str,help="Input toy file")
parser.add_argument("-o", dest="name",default="", type=str,help="Output binned toy file")

args, unknown = parser.parse_known_args()

rt.gROOT.SetBatch(True)
data = rt.TFile.Open(args.data, 'READ')
toy  = rt.TFile.Open(args.toy , 'READ')

ws_data = data.Get('ws_bkg')
data_data = ws_data.obj('data_obs')
data_obs  = ws_data.var('mass_ll')
data_hist = data_data.createHistogram('data_hist', data_obs)

nbins = data_hist.GetNbinsX()
xmin  = data_hist.GetBinLowEdge(1)
xmax  = data_hist.GetXaxis().GetBinUpEdge(nbins)

ws_toy = toy.Get('ws_bkg')
toy_data = ws_toy.obj('data_obs')
toy_obs  = ws_toy.var('mass_ll')
toy_hist = toy_data.createHistogram('toy_hist', toy_obs, rt.RooFit.Binning(nbins, xmin, xmax))

out_data = rt.RooDataHist('data_obs', 'Toy data (binned)', data_obs, toy_hist)

ws_out = rt.RooWorkspace('ws_bkg')
ws_out.Import(out_data)
ws_out.writeToFile(args.name)

