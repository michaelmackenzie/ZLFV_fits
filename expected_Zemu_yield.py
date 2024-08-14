import os
import argparse
import ROOT as rt
import math

def binomial_error(num,den):
  return math.sqrt( num/den*(1-num/den) / den)

def error_mult(vals,errs):
  error=0.
  for ival in range(len(vals)):
    tmp_error=errs[ival]
    for ival2 in range(len(vals)):
      if ival==ival2: continue
      tmp_error*=vals[ival]
    error+=tmp_error*tmp_error
  return math.sqrt(error)

def error_ratio(num,den,dnum,dden):
  return math.sqrt( (1./den*dnum)**2+(num/(den**2)*dden)**2 )


path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/"

files={"2016":path+"Meas_fullAndSF_bdt_v7_signal_mcRun16.root","2017":path+"Meas_fullAndSF_bdt_v7_signal_mcRun17.root","2018":path+"Meas_fullAndSF_bdt_v7_signal_mcRun18.root","Run2":path+"Meas_fullAndSF_bdt_v7_signal_mcRun1*.root"}


bins=[{"xgbmin":"0.3","xgbmax":"0.7"},{"xgbmin":"0.7","xgbmax":"0.9"},{"xgbmin":"0.9","xgbmax":"1.01"}]
ndens={"2016":180400.,"2017":187000.,"2018":194000.,"Run2":561400.}
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":146.9}
sf="Electron_RecoID_wt*Muon_ID_wt*Electron_ID_wt*Muon_IsoID_wt*Electron_IsoID_wt*PU_wt*PtZ_wt*PtSignal_wt*Trg_wt"
BF=2.62e-7
xsec=6077220.0/(3*0.0336)


rt.gROOT.SetBatch(True)

for year in ["2016","2017","2018","Run2"]:
  cc=rt.TChain("mytreefit")
  cc.Add(files[year])
  print year,"total MC:",cc.GetEntries()
  
  for bn in bins:
    print " - bin",bn
    htemp = rt.TH1F("signal_"+bn['xgbmin']+"_"+bn['xgbmax'],"",25,70,110)
    cc.Draw("mass_ll>>"+"signal_"+bn['xgbmin']+"_"+bn['xgbmax'],sf+"*(70<mass_ll && mass_ll<110 && "+bn['xgbmin']+"<xgb && xgb<"+bn['xgbmax']+")")
    nnum = htemp.Integral()
    nsgn = lumis[year]*xsec*BF*nnum/ndens[year]
    print "   -- a*e ",nnum/ndens[year],"+/-",binomial_error(nnum,ndens[year])
    print "   -- yield ",nsgn,"+/-",lumis[year]*xsec*BF*binomial_error(nnum,ndens[year])

    
  
