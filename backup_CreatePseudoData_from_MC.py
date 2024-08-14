import ROOT as rt
from array import array


samples = ["/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSFAndGenDecay_bdt_v7_bkg_dy_mcRun1*.root",\
"/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndGenWTonly_bdt_v7_bkg_ww_mcRun1*.root",\
"/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndGenWTonly_bdt_v7_bkg_ttbar_mcRun1*.root",\
"/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root"]


norm_samples = [ 
    ["/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2018MCDATA/MC/DY50_2l_MC18/LFVAnalysis_DY50-amc_2018.root",\
     "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2017MCDATA/MC/DY50_2l_MC17/LFVAnalysis_DY50-ext_2017.root",\
     "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2016MCDATA/MC/DY50_2l_MC16/LFVAnalysis_DY50-amc_2016.root" ],\
   
   ["/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2016MCDATA/MC/WW2lnu_MC16/LFVAnalysis_WW_2016.root",\
    "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2017MCDATA/MC/WW2lnu_MC17/LFVAnalysis_WW_2017.root",\
    "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2018MCDATA/MC/WW2lnu_MC18/LFVAnalysis_WW_2018.root"],\

  ["/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2016MCDATA/MC/TTbarLeptonic_MC16/LFVAnalysis_ttbarlnu_2016.root",\
   "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2017MCDATA/MC/TTbarLeptonic_MC17/LFVAnalysis_ttbarlnu_2017.root",\
   "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/2018MCDATA/MC/TTbarLeptonic_MC18/LFVAnalysis_ttbarlnu_2018.root"],\

  ["/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/LFVAnalysis_ZEMu-v2_2016.root",\
   "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/LFVAnalysis_ZEMu-v2_2017.root",\
   "/eos/cms/store/group/phys_smp/ZLFV/lfvanalysis_rootfiles/LFVAnalysis_ZEMu-v2_2018.root"]

]
#data_sample_for_weight = "/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root"
Br_xsecs=[6077220.0,12178.0,88310.0,15.79] #6077220.0/(3*0.0336)*2.62e-7]
sample_names=["dy","ww","ttbar","signal"]
denominator_fraction=[1,0.5,0.3,1]
dy_index=0
signal_index=3
signal_r=2



#datatree = rt.TTree("mytreefit","mytreefit")
#datatree.Add(data_sample_for_weight)



outtree = rt.TTree("mytreefit","mytreefit")

bdt = array('f',[0])
mass = array('f',[0])
sources = [ array('f',[0]) for i in range(len(sample_names)) ]
Gen_wt = array('f',[0])
IsZmm = array('f',[0])
IsZee = array('f',[0])
IsZtt = array('f',[0])
JetPUIDWeight = array('f',[0])
SignGen_wt = array('f',[0])
MixZ_wt = array('f',[0])
SFbtag_wt = array('f',[0])
Muon_dxydz_wt = array('f',[0])
Electron_dxydz_wt = array('f',[0])
Trg_wt = array('f',[0])
PtSignal_wt = array('f',[0])
PtZ_wt = array('f',[0])
Prefire_wt = array('f',[0])
PU_wt = array('f',[0])
Total_wt = array('f',[0])
Muon_RecoID_wt = array('f',[0])
Electron_RecoID_wt = array('f',[0])
Muon_ID_wt = array('f',[0])
Electron_ID_wt = array('f',[0])
Muon_IsoID_wt = array('f',[0])
Electron_IsoID_wt = array('f',[0])
Flag_met = array('f',[0])
Flag_muon = array('f',[0])
Flag_electron = array('f',[0])

outtree.Branch("xgb",bdt,"xgb/F")
outtree.Branch("mass_ll",mass,"mass_ll/F")
for i in range(len(sources)):
  outtree.Branch(sample_names[i],sources[i],sample_names[i]+"/F")
outtree.Branch("IsZmm",IsZmm,"IsZmm/F")
outtree.Branch("IsZee",IsZee,"IsZee/F")
outtree.Branch("IsZtt",IsZtt,"IsZtt/F")
outtree.Branch("Gen_wt",Gen_wt,"Gen_wt/F")
outtree.Branch("SignGen_wt",SignGen_wt,"SignGen_wt/F")


outtree.Branch("JetPUIDWeight",JetPUIDWeight,"JetPUIDWeight/F")
outtree.Branch("MixZ_wt",MixZ_wt,"MixZ_wt/F")
outtree.Branch("SFbtag_wt",SFbtag_wt,"SFbtag_wt/F")
outtree.Branch("Trg_wt",Trg_wt,"Trg_wt/F")
outtree.Branch("PtSignal_wt",PtSignal_wt,"PtSignal_wt/F")
outtree.Branch("PtZ_wt",PtZ_wt,"PtZ_wt/F")
outtree.Branch("Prefire_wt",Prefire_wt,"Prefire_wt/F")
outtree.Branch("PU_wt",PU_wt,"PU_wt/F")
outtree.Branch("Total_wt",Total_wt,"Total_wt/F")
outtree.Branch("Flag_met",Flag_met,"Flag_met/F")
outtree.Branch("Flag_muon",Flag_muon,"Flag_muon/F")
outtree.Branch("Flag_electron",Flag_electron,"Flag_electron/F")

nums=[]
dens=[]

rt.gROOT.SetBatch(True)
for ismp in range(len(samples)):
   intree_num = rt.TChain("mytreefit")
   intree_num.Add(samples[ismp])
   intree_den = rt.TChain("Norm")
   intree_den.Add(norm_samples[ismp][0])
   intree_den.Add(norm_samples[ismp][1])
   intree_den.Add(norm_samples[ismp][2])
   den=0
   for evt in intree_den:
      den+=(evt.NEvents-evt.NNegative)*denominator_fraction[ismp]
   htemp = rt.TH1F("hnum_"+str(ismp),"",1,0,2)
   intree_num.Draw("1>>hnum_"+str(ismp),"SignGen_wt")
   num = htemp.GetBinContent(1)
   nums.append(num)
   dens.append(den)

       
dy_evts=1
for ismp in range(len(samples)):
   intree = rt.TChain("mytreefit")
   intree.Add(samples[ismp])
   print "sample",sample_names[ismp],"in root evts",intree.GetEntries()
   print "weighted num",nums[ismp],"den",dens[ismp]
   if ismp==dy_index:
      dy_evts=intree.GetEntries()
   evt_max=nums[ismp]/dens[ismp]*Br_xsecs[ismp]/(nums[0]/dens[0]*Br_xsecs[0])*dy_evts
   if signal_index==ismp:
     evt_max=intree.GetEntries()
     print "will run on",evt_max,"normalized to",signal_r*(nums[ismp]/dens[ismp]*Br_xsecs[ismp]/(nums[0]/dens[0]*Br_xsecs[0]))*nums[0]
   else:
     print "will run on",evt_max
   ievt=0
   for evt in intree:
      if ievt>int(evt_max): break;
      for idef in range(len(samples)): 
          sources[idef][0]=0
      bdt[0] = evt.xgb
      sources[ismp][0]=1
      mass[0]=evt.mass_ll
      if signal_index==ismp:
        SignGen_wt[0] = signal_r*(nums[ismp]/dens[ismp]*Br_xsecs[ismp]/(nums[0]/dens[0]*Br_xsecs[0]))*nums[0]/evt_max
      else:
        SignGen_wt[0] = evt.SignGen_wt
      Gen_wt[0] = evt.Gen_wt
      Flag_met[0]=evt.Flag_met
      Flag_muon[0]=evt.Flag_muon
      Flag_electron[0]=evt.Flag_electron
      if ismp==dy_index:
         IsZmm[0]=evt.IsGen_Zmm
         IsZee[0]=evt.IsGen_Zee
         IsZtt[0]=evt.IsGen_Ztt
      else:
        IsZmm[0], IsZee[0], IsZtt[0] = 0, 0, 0
#      JetPUIDWeight[0]=evt.JetPUIDWeight
      '''
      MixZ_wt[0]=evt.MixZ_wt
      SFbtag_wt[0]=evt.SFbtag_wt
      Trg_wt[0]=evt.Trg_wt
      PtSignal_wt[0]=evt.PtSignal_wt
      PtZ_wt[0]=evt.PtZ_wt
      Prefire_wt[0]=evt.Prefire_wt
      PU_wt[0]=evt.PU_wt'''
      #Total_wt[0]=JetPUIDWeight[0]*SignGen_wt[0]*MixZ_wt[0]*SFbtag_wt[0]*Trg_wt[0]*PtSignal_wt[0]*PtZ_wt[0]*Prefire_wt[0]*Gen_wt[0]/abs(Gen_wt[0])*PU_wt[0]
      ievt+=1
      outtree.Fill()
         
fout = rt.TFile("pseudo_data_from_MC_r"+str(signal_r)+".root","RECREATE")
outtree.Write()
