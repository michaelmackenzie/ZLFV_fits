import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TGraphAsymmErrors
from array import array
import sys
sys.path.insert(1, '/afs/cern.ch/work/g/gkaratha/private/Analysis/DispJets/Analyzer/CMSSW_10_2_16_UL/src/PhysicsTools/NanoAODTools/plotter/ZMuE_plotting_and_cfg/BDT/')
from tdrstyle import *
from cms_lumi import *
# from Plot_helper import *
ROOT.gROOT.SetBatch(True)


BF=2.62e-7
categories=array('f',[1,2,3,4])
labels=["BDT bin1","BDT bin2","BDT bin3","All BDT bins"]
name="limits"

# --------------------- v01 results ---------------------------------
#expected.
# r_exp_low_2sigma=[1.2638,0.6188,0.6236,0.4255]
# r_exp_low_1sigma=[1.7639, 0.8093,0.8286,0.5568]
# r_exp_central=[2.7188, 1.0703,1.1484,0.7617]
# r_exp_high_1sigma=[3.6832, 1.7785,2.2013,1.0502]
# r_exp_high_2sigma=[4.9188,2.8280,3.1392,1.4706]
# #observed
# #r_obs_central=[2.2886, 1.2077,1.3957,0.7938] #data
# r_obs_central=[2.7188, 1.0703,1.1484,0.7617] #asimov

# --------------------- v02 results ---------------------------------
#expected
BF = 1.e-7
name = "limits_v02"
r_exp_low_2sigma =[5.0476 , 1.6259, 1.6326, 1.0971]
r_exp_low_1sigma =[6.2942 , 2.1153, 2.1655, 1.4465]
r_exp_central    =[8.2469 , 2.8867, 3.0055, 2.0398]
r_exp_high_1sigma=[10.7739, 4.1165, 4.2020, 2.8418]
r_exp_high_2sigma=[13.7387, 5.4799, 5.6508, 3.8478]
#observed
r_obs_central    =[5.7365 , 3.1833, 3.0038, 1.9198]

## upd
#r_exp_low_2sigma = [1.2602, 0.6130, 0.6299, 0.4244]
#r_exp_low_1sigma = [1.9289, 0.8030, 0.8370, 0.5554]
#r_exp_central = [2.7109, 1.0898, 1.1602, 0.7598]
#r_exp_high_1sigma = [3.6294, 1.9847, 2.2423, 1.0475]
#r_exp_high_2sigma = [4.8911, 2.7844, 3.1562, 1.7036]
#observed
#r_obs_central = [2.1647, 1.2085, 1.3191, 0.7336] #data
#r_obs_central=[2.7109, 1.0898, 1.1484, 0.7598] #asimov

###############################################################################



exp_central=array('f',[])
exp_high_sigma1=array('f',[])
exp_high_sigma2=array('f',[])
exp_low_sigma1=array('f',[])
exp_low_sigma2=array('f',[])

obs_central=array('f',[])


for idx in range(len(r_exp_central)):
  exp_central.append(BF*r_exp_central[idx])
  exp_high_sigma1.append(BF*(r_exp_high_1sigma[idx]-r_exp_central[idx]))
  exp_high_sigma2.append(BF*(r_exp_high_2sigma[idx]-r_exp_central[idx]))
  exp_low_sigma1.append(BF*(r_exp_central[idx]-r_exp_low_1sigma[idx]))
  exp_low_sigma2.append(BF*(r_exp_central[idx]-r_exp_low_2sigma[idx]))
  obs_central.append(BF*r_obs_central[idx])
  print(labels[idx],"exp",exp_central[-1],"obs",obs_central[-1])



def plotUpperLimits(categories,exp_central,exp_low_1sigma,exp_high_1sigma,exp_low_2sigma,exp_high_2sigma,obs_central,name):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
     
    N = len(categories)
 
    up2s = []
    
    x_error = array("f")
    y_error = array("f")
    for i in range(N):
      up2s.append(exp_central[i]+exp_high_1sigma[i])
      x_error.append(0.5)
      y_error.append(0.01*exp_central[i])

    exp_2sigma = TGraphAsymmErrors(N,categories,exp_central,x_error,x_error,exp_low_sigma2,exp_high_sigma2)
    exp_1sigma = TGraphAsymmErrors(N,categories,exp_central,x_error,x_error,exp_low_sigma1,exp_high_sigma1)    
    exp_median = TGraphAsymmErrors(N,categories,exp_central,x_error,x_error,y_error,y_error)
    obs_median = TGraphAsymmErrors(N,categories,obs_central,x_error,x_error,y_error,y_error)
    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()
    frame = c.DrawFrame(-0.1,0.001, max(categories)+0.7, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on BF(Z^{0} #rightarrow e^{#pm} #mu^{#mp})")
    frame.GetXaxis().SetTitle("")
    frame.SetMinimum(0)
    frame.SetMaximum(max(up2s)*1.75)
#    for cat,label in zip(categories,labels):
#      print int(cat),label
#      frame.GetXaxis().SetBinLabel(int(cat),label)
    frame.GetXaxis().ChangeLabel(1,-1,-1,-1,-1,-1," ")
    for cat,label in zip(categories,labels):
      frame.GetXaxis().ChangeLabel(int(1+cat),-1,-1,-1,-1,-1,label) 
    
    exp_2sigma.SetFillColor(ROOT.kOrange)
    exp_2sigma.SetLineColor(ROOT.kOrange)
    exp_2sigma.SetMarkerSize(0)
    exp_2sigma.SetFillStyle(1001)
    exp_2sigma.Draw('F2')
 
    exp_1sigma.SetFillColor(ROOT.kGreen+1)
    exp_1sigma.SetLineColor(ROOT.kGreen+1)
    exp_1sigma.SetMarkerSize(0)
    exp_1sigma.SetFillStyle(1001)
    exp_1sigma.Draw('F2')
 
    exp_median.SetLineColor(1)
    exp_median.SetLineWidth(2)
    exp_median.SetMarkerStyle(0)
    exp_median.SetMarkerSize(0)
    exp_median.SetLineStyle(2)
    exp_median.Draw('Psame')

    obs_median.SetLineColor(0)
    obs_median.SetLineWidth(0)
    obs_median.SetMarkerStyle(8)
    obs_median.SetMarkerSize(1)
    obs_median.Draw('Psame')
 
    CMS_lumi(c,14,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')
#    frame.GetXaxis().SetTitle("BDT categories") 
    x1 = 0.55
    x2 = x1 + 0.24
    y2 = 0.76
    y1 = 0.60
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(obs_median, "Asymptotic CL_{s} observed",'P')
    legend.AddEntry(exp_median, "Asymptotic CL_{s} expected",'L')
    legend.AddEntry(exp_1sigma, "#pm 1 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 1 std. deviation",'f')
    legend.AddEntry(exp_2sigma,"#pm 2 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 2 std. deviation",'f')
    legend.Draw()
 
    
    c.SaveAs("lp_"+name+".png")
    c.Close()
 
 
# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step

if __name__ == '__main__':
   plotUpperLimits(categories,exp_central,exp_low_sigma1,exp_high_sigma1,exp_low_sigma2,exp_high_sigma2,obs_central,name)
