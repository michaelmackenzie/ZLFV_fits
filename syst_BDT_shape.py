import os
import argparse
import ROOT as rt
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')


run="Run2"
systematics=["MU_SCALEDOWN","MU_SCALEUP","ELE_SCALEUP","ELE_SCALEDOWN","MET_JER","MET_JES"]
mass_systematic={"MU_SCALEDOWN":"mass_ll_Muon_scale_down","MU_SCALEUP":"mass_ll_Muon_scale_up","ELE_SCALEUP":"mass_ll_Electron_scale_up","ELE_SCALEDOWN":"mass_ll_Electron_scale_down","MET_JER":"mass_ll","MET_JES":"mass_ll"}
path="../BDT/Systematics_v2_no_met_significance/"
nominals=[122.53, 219.88, 79.92]

bins=[{"xgbmin":"0.3","xgbmax":"0.7"},{"xgbmin":"0.7","xgbmax":"0.9"},{"xgbmin":"0.9","xgbmax":"1.01"}]
ndens={"2016":180400.,"2017":187000.,"2018":194000.,"Run2":561400.}
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}
sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)"
BF=2.62e-7
xsec=6077220.0/(3*0.0336)

yields={syst:[] for syst in systematics}



for syst in systematics:   
  print "Systematic:",syst
  files={"2016":path+"/Meas_fullAndSF_syst_"+syst+"_bdt_v7_signal_mcRun16.root","2017":path+"/Meas_fullAndSF_syst_"+syst+"_bdt_v7_signal_mcRun17.root","2018":path+"/Meas_fullAndSF_syst_"+syst+"_bdt_v7_signal_mcRun18.root"}
  chains={}
  for year in files.keys():
    chains[year]=rt.TChain("mytreefit")
    chains[year].Add(files[year])
 
  if run =="Run2":
    for bn in bins:
      nsgn=0.
      for year in files.keys():  
        htemp = rt.TH1F("signal_"+syst+"_"+bn['xgbmin']+"_"+bn['xgbmax']+"_"+year,"",25,70,110)
        chains[year].Draw(mass_systematic[syst]+">>"+"signal_"+syst+"_"+bn['xgbmin']+"_"+bn['xgbmax']+"_"+year,sf+"*(70<"+mass_systematic[syst]+" && "+mass_systematic[syst]+"<110 && "+bn['xgbmin']+"<xgb && xgb<"+bn['xgbmax']+" && Flag_met && Flag_muon && Flag_electron)")
        nnum = htemp.Integral()
        nsgn += lumis[year]*xsec*BF*nnum/ndens[year]
      print "  - ",nsgn 
      yields[syst].append(nsgn)
  else:
    nden=ndens[run]
    cc=rt.TChain("mytreefit")
    cc.Add(files[run])
    for bn in bins:
      print "  bin",bn
      htemp = rt.TH1F("signal_"+syst+"_"+bn['xgbmin']+"_"+bn['xgbmax'],"",25,70,110)
      cc.Draw(mass_systematic[syst]+">>"+"signal_"+syst+"_"+bn['xgbmin']+"_"+bn['xgbmax'],sf+"*(70<"+mass_systematic[syst]+" && "+mass_systematic[syst]+"<110 && "+bn['xgbmin']+"<xgb && xgb<"+bn['xgbmax']+" && Flag_met && Flag_muon && Flag_electron)")
      nnum = htemp.Integral()
      nsgn = lumis[run]*xsec*BF*nnum/ndens[run]
      print "  - ",nsgn
      yields[syst].append(nsgn)

for ibin in range(len(bins)):
  print "bin",ibin,bins[ibin],"systematics:"
  for syst in systematics:
    print "  - ",syst,"ratio",yields[syst][ibin]/nominals[ibin]," or ",abs(yields[syst][ibin]-nominals[ibin])/nominals[ibin]*100,"%"
   
    
  
