import ROOT as rt
import math
#rt.gInterpreter.Declare('#include "SFBDT_weight.h"')
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')


rt.gROOT.SetBatch(True)

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


path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/NoMET_significance/"
years =['2016','2017','2018']
files={"2016":"Meas_fullAndSF_bdt_v7_signal_mcRun16_no_met_significance.root","2017":"Meas_fullAndSF_bdt_v7_signal_mcRun17_no_met_significance.root","2018":"Meas_fullAndSF_bdt_v7_signal_mcRun18_no_met_significance.root"}
ndens={"2016":180400.,"2017":187000.,"2018":194000.}
lumis={"2016":36.33,"2017":41.48,"2018":59.83}
total_lumi= 137.6
xsecBF = 6077220.0/(3*0.0336)*2.62e-7
systematics=["nominal","recoEle","idMu","idEle","isoMu","isoEle","dxyMu","dxyEle","pu","trg","ptz","prefire","btag","mix","pdfsys"]
Debug=False



############################# systematics detailed cfg ########################
bins={"bin1":"(Flag_met && Flag_muon && Flag_electron && 70<mass_ll && mass_ll<110 && 0.3<xgb && xgb<0.7)",\
      "bin2":"( Flag_met && Flag_muon && Flag_electron && 70<mass_ll && mass_ll<110 && 0.7<xgb && xgb<0.9)",\
      "bin3":"( Flag_met && Flag_muon && Flag_electron && 70<mass_ll && mass_ll<110 && 0.9<xgb )"}

variations={"nominal":[""],"recoEle":["UP","DOWN"],"idMu":["UP","DOWN"],"idEle":["UP","DOWN"],"isoMu":["UP","DOWN"],"isoEle":["UP","DOWN"],"ptz":[""],"dxyMu":['UP','DOWN'],"dxyEle":['UP','DOWN'],"pu":["UP","DOWN"],"trg":["UP","DOWN"],"prefire":["UP","DOWN"],"btag":["UP","DOWN"],"mix":["UP","DOWN"],"pdfsys":[""]}

sfs={ 'nominal':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'recoEleUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_up*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'recoEleDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_down*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'idMuUP':'Muon_RecoID_wt*Muon_ID_up*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'idMuDOWN':'Muon_RecoID_wt*Muon_ID_down*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'idEleUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_up*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'idEleDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_down*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'isoMuUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_up*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'isoMuDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_down*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'isoEleUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_up*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'isoEleDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_down*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'dxyMuUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_up*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'dxyMuDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_down*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'dxyEleUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_up*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'dxyEleDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_down*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'puUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_up*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'puDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_down*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'ptz':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_sys*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'trgUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_up*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'trgDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_down*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'prefireUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_up*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
      'prefireDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_down*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
       'btagUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_up*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
       'btagDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_down*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
       'mixUP':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_up*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
       'mixDOWN':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_down*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*',\
       'pdfsys':'Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)*PDF_sys*'
}
     

###############################################################################
################################# code ########################################
###############################################################################


yields={i:{j:{k:0 for k in systematics} for j in years} for i in bins.keys()}

stat_error={i:{k:0 for k in years} for i in bins.keys()}

chains={}
for year in years:
  chains[year]=rt.TChain("mytreefit")
  chains[year].Add(path+"/"+files[year])



for bn in bins.keys():
  if Debug: print "\nbin:",bn
  for syst in systematics:
    if Debug: print " -Systematic:",syst
    for year in years:
      if Debug: print "  --year:",year
      variation_yields=[]
      for var in variations[syst]:
        if Debug: print "  ---Variation:",var
        sf_temp=sfs[syst+var]
        hnum = rt.TH1F("hnum_"+syst+var+year+"_"+bn,"",1,0,2)
        chains[year].Draw("1>>hnum_"+syst+var+year+"_"+bn,sf_temp+bins[bn])
        num=hnum.Integral()
        if Debug:
           print "    ---->>> ae=",num/ndens[year],"+/-",binomial_error(num,ndens[year]),"Nsignal=",lumis[year]*xsecBF*num/ndens[year],"+/-",lumis[year]*xsecBF*binomial_error(num,ndens[year])
        if syst=="nominal":
           stat_error[bn][year]= lumis[year]*xsecBF*binomial_error(num,ndens[year])
        variation_yields.append(lumis[year]*xsecBF*num/ndens[year])

      if len(variation_yields)==1:
         yields[bn][year][syst]=variation_yields[0]
      else:
         yields[bn][year][syst]=variation_yields[0] if abs(variation_yields[0]-yields[bn][year]['nominal']) > abs(variation_yields[1]-yields[bn][year]['nominal']) else variation_yields[1]

run2_yields ={i:{k:0 for k in systematics} for i in bins.keys()}
run2_stat_error= {i:0 for i in bins.keys()}


for k in systematics:
  for i in bins.keys():
    run2_yields[i][k] =  yields[i]['2016'][k]+yields[i]['2017'][k]+yields[i]['2018'][k]

for i in bins.keys():
  run2_stat_error[i] = math.sqrt(stat_error[i]['2016']**2+stat_error[i]['2017']**2+stat_error[i]['2018']**2)

print "results"
for bn in bins.keys():
  print " - bin",bn  
  for syst in systematics:
    print "  -- "+syst+"=",run2_yields[bn][syst],syst+"/nominal=",(run2_yields[bn][syst])/run2_yields[bn]['nominal']," or ",abs(run2_yields[bn][syst]-run2_yields[bn]['nominal'])/run2_yields[bn]['nominal']*100,"%"
  print "  -- statistic error",run2_stat_error[bn],"ratio",1+run2_stat_error[bn]/run2_yields[bn]['nominal'] 
