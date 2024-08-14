import ROOT as rt
import math
#rt.gInterpreter.Declare('#include "SFBDT_weight.h"')

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


path_emu="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun*.root"
nden_emu=561400.
path_mumu="Meas_bdt_v5_Zmumu_mcRun18_v3_part*.root"
nden_mumu=64401410.
path_ee="Meas_bdt_v5_Zee_mcRun18_v3_part*.root"
nden_ee=64400859.  

cuts={"mass":"(70<mass_ll && mass_ll<110)","bdt":"(xgb>0.3 && xgb<0.7)","all":"(70<mass_ll && mass_ll<110 && xgb>0.3 && xgb<0.7 )"} #all MUST be always here

sf="Muon_RecoID_wt*Electron_RecoID_wt*Muon_ID_wt*Electron_ID_wt*Muon_IsoID_wt*Electron_IsoID_wt*PtZ_wt*PU_wt*Prefire_wt*Trg_wt*"
sf_only_emu="PtSignal_wt*"
sf_bdt="" #no "*" in end; puts directly the f state deactivate with ""
list_eff_samples=["emu"] # max: "mumu","ee","emu"
Calculate_ratio=False

###############################################################################
###############################################################################
###############################################################################
paths={"emu":path_emu,"ee":path_ee,"mumu":path_mumu}
ndens={"emu":nden_emu,"ee":nden_ee,"mumu":nden_mumu}
ae_value={"emu":-1.,"ee":-1.,"mumu":-1.}
ae_error={"emu":-1.,"ee":-1.,"mumu":-1.}

for smpl in list_eff_samples:
  cc=rt.TChain("mytreefit")
  cc.Add(paths[smpl])
  print "\nprinting ae for Z->",smpl,":"
  den_total=ndens[smpl]
  hden = rt.TH1F("hden_"+smpl,"",1,0,2)
  sf_temp=sf
  if "emu"==smpl: sf_temp+=sf_only_emu
  if sf_bdt!="": sf_temp+=sf_bdt+smpl+"(xgb)*"
  cc.Draw("1>>hden_"+smpl,sf_temp+"(1)")
  den_reco=hden.Integral()
  print "- ae(reco)",den_reco/den_total,"+/-",binomial_error(den_reco,den_total)
  print "- individual eff for cut:"
  num_all=0.
  for cut in cuts.keys(): 
    hnum = rt.TH1F("hnum_"+smpl+"_"+cut,"",1,0,2)
    cc.Draw("1>>hnum_"+smpl+"_"+cut,sf_temp+cuts[cut])
    num=hnum.Integral()
    if "all"==cut:
      num_all=num
    else:
      print "  -- ae(",cut,")",num/den_reco,"+/-",binomial_error(num,den_reco)

  print "- ae(sel-ALL cuts)",num_all/den_reco,"+/-",binomial_error(num_all,den_reco),num_all,den_reco
  print "-> a*e(total)=",num_all/den_total,"+/-",binomial_error(num_all,den_total)
  ae_value[smpl]=num_all/den_total
  ae_error[smpl]=binomial_error(num_all,den_total)


if Calculate_ratio and all(ae_value[x]>-1 for x in ae_value.keys()):
   num_val = ae_value["mumu"]*ae_value["ee"]
   sqt_num_val = math.sqrt(num_val)
   sqt_num_error =math.sqrt( (math.sqrt(ae_value["mumu"]/ae_value["ee"])*ae_error["ee"]/2.0 )**2 + (math.sqrt(ae_value["ee"]/ae_value["mumu"])*ae_error["mumu"]/2.0)**2 )
   ratio = sqt_num_val/ae_value["emu"]
   ratio_error = error_ratio(sqt_num_val, ae_value["emu"], sqt_num_error, ae_error["emu"])
   print "\n----------------------------"
   print "ratio",ratio,"+/-",ratio_error
   print "----------------------------\n"   

