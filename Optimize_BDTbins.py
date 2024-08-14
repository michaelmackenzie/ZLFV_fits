import 	ROOT as rt
from ROOT import RooFit
import os
rt.gInterpreter.Declare('#include "SFBDT_weight.h"')


def Getting_MC_exp(path,nden,sf,bdt_min,bdt_max,BFxLxSigma):
   cc=rt.TChain("mytreefit")
   cc.Add(path)
   htemp = rt.TH1F("htemp_"+str(bdt_min)+"_"+str(bdt_max),"",1,0,2)   
   cc.Draw("1>>htemp_"+str(bdt_min)+"_"+str(bdt_max),sf+"*(85<mass_ll && mass_ll<95 && "+str(bdt_min)+"<xgb && xgb<"+str(bdt_max)+")")
   nnum = htemp.Integral()
   return BFxLxSigma*nnum/nden

    
def Getting_Bks_inSR(path,name,bdt_min,bdt_max):
   msgservice = rt.RooMsgService.instance()
   msgservice.setGlobalKillBelow(RooFit.FATAL)
   cc=rt.TChain("mytreefit")
   cc.Add(path)
   cc_cut = cc.CopyTree("70<mass_ll && mass_ll<110 && "+str(bdt_min)+"<xgb && xgb<"+str(bdt_max),"")
   dilep_mass = rt.RooRealVar("mass_ll","m(e,#mu)", 110, 70, 110)
   dataset = rt.RooDataSet("data","data",cc_cut,rt.RooArgSet(dilep_mass))
   dilep_mass.setRange("window",85,95)
   bkg_exp_alpha = rt.RooRealVar("bkg_exp_alpha", "bkg_exp_alpha", 1., -10., 10.)
   bkg_exp = rt.RooExponential("bkg_exp", "bkg_exp", dilep_mass, bkg_exp_alpha)
   nbkg = rt.RooRealVar("nbkg","",1000, 0, 1000000)
   ebkg_exp = rt.RooExtendPdf("ebkg_exp","ebkg_exp",bkg_exp,nbkg)
   bkg_result = ebkg_exp.fitTo(dataset,RooFit.Extended(True),RooFit.Save(), RooFit.PrintLevel(-1))
   bkg_result.Print()
   xframe = dilep_mass.frame()
   dataset.plotOn(xframe, RooFit.Binning(50)) 
   ebkg_exp.plotOn(xframe)
   cnv = rt.TCanvas("opt_test_"+str(bdt_min)+"_"+str(bdt_max), "CMS Preliminary", 700, 700)
   xframe.Draw()
   cnv.SaveAs("opt_bkg_"+name+"_min_"+str(bdt_min)+"_max_"+str(bdt_max)+".png")
   return nbkg.getVal()*ebkg_exp.createIntegral(rt.RooArgSet(dilep_mass),rt.RooArgSet(dilep_mass),"window").getVal()



if __name__ == "__main__":
  rt.gROOT.SetBatch(True)
  #getting the signal
  path_emu="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v5_tuples/BDT_outputs_v5/Meas_fullAndSF_bdt_v5_signal_mcRun18.root"
  nden_emu=0.2*194000.
  sf_emu="recol1_weight*recol2_weight*idl1_weight*idl2_weight*isol1_weight*isol2_weight*pu_weight*zpt_weight*trigger_weight*zsgn_weight*SFBDT_weight_Zemu(xgb)"
  BFxLxSigma=64*6435000.0/(3*0.0336)*2.62e-7
  path_data_emu="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v5_tuples/BDT_outputs_v5/Meas_full_bdt_v5_data_emu_Run18.root"
  run_1bin=False
  run_2bin=False
  run_3bin=True
  run_4bin=False
  extra_name="_0p5" #start with "_" to be better visually
  start_val=0.5
  step=0.1




  start_bin=int(start_val/step)

  if run_1bin:
    os.system("rm result_combine_1bin"+extra_name+".txt")
    for icut in range(start_bin,10):
       # 1 Bin analysis
       datacard_1bin="imax    1 number of bins\njmax    1 number of processes minus 1\nkmax    * number of nuisance parameters\n"
       datacard_1bin+="----------------------------\n---------------------------\n"
       datacard_1bin+="bin signal_region\n"
       bdt_cut = icut*step
       sgn= Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut,star_val,BFxLxSigma) 
       bkg= Getting_Bks_inSR(path_data_emu,"bin1",bdt_cut,1.0)
       datacard_1bin+="observation "+str(bkg)+"\n"
       datacard_1bin+="bin signal_region signal_region\n"
       datacard_1bin+="process hbackground hsignal\n"
       datacard_1bin+="process 1 0\n"
       datacard_1bin+="rate "+str(bkg)+" "+str(sgn)
       with open("combine_bin1_cut_"+str(bdt_cut)+extra_name+".txt",'w') as txt:
          txt.write(datacard_1bin)
          txt.close()
       os.system('echo "bdt_cut '+str(bdt_cut)+extra_name+'">>result_combine_1bin'+extra_name+'.txt')
       os.system("./run_combine.sh combine_bin1_cut_"+str(bdt_cut)+extra_name+".txt combine_bin1_cut_"+str(bdt_cut)+extra_name+" >> result_combine_1bin"+extra_name+".txt")

  # 2 Bin analysis      
  if run_2bin:
    os.system("rm result_combine_2bin"+extra_name+".txt")
    for icut in range(start_bin,10):
       # 2 Bin analysis
       datacard_2bin="imax    2 number of bins\njmax    * number of processes minus 1\nkmax    * number of nuisance parameters\n"
       datacard_2bin+="----------------------------\n---------------------------\n"
       datacard_2bin+="bin signal_region1 signal_region2\n"
       bdt_cut = icut*step
       sgn1 = Getting_MC_exp(path_emu,nden_emu,sf_emu,start_val,bdt_cut,BFxLxSigma) 
       bkg1= Getting_Bks_inSR(path_data_emu,"bin2a",start_val,bdt_cut)
       sgn2 = Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut,1.0,BFxLxSigma)    
       bkg2 = Getting_Bks_inSR(path_data_emu,"bin2b",bdt_cut,1.0)
       datacard_2bin+="observation "+str(bkg1)+" "+str(bkg2)+"\n"
       datacard_2bin+="bin signal_region1 signal_region1 signal_region2 signal_region2\n"
       datacard_2bin+="process background signal background signal\n"
       datacard_2bin+="process 1 0 1 0\n"
       datacard_2bin+="rate "+str(bkg1)+" "+str(sgn1)+"  "+str(bkg2)+" "+str(sgn2)
       with open("combine_bin2_cut_"+str(bdt_cut)+extra_name+".txt",'w') as txt:
          txt.write(datacard_2bin)
          txt.close()
       os.system('echo bdt_cut '+str(bdt_cut)+extra_name+'>>result_combine_2bin'+extra_name+'.txt')
       os.system("./run_combine.sh combine_bin2_cut_"+str(bdt_cut)+extra_name+".txt combine_bin2_cut_"+str(bdt_cut)+extra_name+" >> result_combine_2bin"+extra_name+".txt")
     
  # 3 Bin analysis      
  if run_3bin:
    os.system("rm result_combine_3bin"+extra_name+".txt")
    for icut1 in range(start_bin,9):
       for icut2 in range(icut1+1,10):
          # 3 Bin analysis
          datacard_3bin="imax    3 number of bins\njmax    * number of processes minus 1\nkmax    * number of nuisance parameters\n"
          datacard_3bin+="----------------------------\n---------------------------\n"
          datacard_3bin+="bin signal_region1 signal_region2 signal_region3\n"
          bdt_cut1 = icut1*step
          bdt_cut2 = icut2*step
          sgn1 = Getting_MC_exp(path_emu,nden_emu,sf_emu,start_val,bdt_cut1,BFxLxSigma) 
          bkg1= Getting_Bks_inSR(path_data_emu,"bdt3a",start_val,bdt_cut1)
          sgn2 = Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut1,bdt_cut2,BFxLxSigma)    
          bkg2 = Getting_Bks_inSR(path_data_emu,"bdt3b",bdt_cut1,bdt_cut2)
          sgn3 = Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut2,1.0,BFxLxSigma)
          bkg3 = Getting_Bks_inSR(path_data_emu,"bdt3c",bdt_cut2,1.0)
          datacard_3bin+="observation "+str(bkg1)+" "+str(bkg2)+" "+str(bkg3)+"\n"
          datacard_3bin+="bin signal_region1 signal_region1 signal_region2 signal_region2 signal_region3 signal_region3\n"
          datacard_3bin+="process background signal background signal background signal\n"
          datacard_3bin+="process 1 0 1 0 1 0\n"
          datacard_3bin+="rate "+str(bkg1)+" "+str(sgn1)+"  "+str(bkg2)+" "+str(sgn2)+" "+str(bkg3)+" "+str(sgn3)
          with open("combine_bin3_cut1_"+str(bdt_cut1)+"_cut2_"+str(bdt_cut2)+extra_name+".txt",'w') as txt:
             txt.write(datacard_3bin)
             txt.close()
          os.system('echo bdt_cut '+str(bdt_cut1)+' bdt_cut2 '+str(bdt_cut2)+'>>result_combine_3bin'+extra_name+'.txt')
          os.system("./run_combine.sh combine_bin3_cut1_"+str(bdt_cut1)+"_cut2_"+str(bdt_cut2)+extra_name+".txt combine_bin3_cut1_"+str(bdt_cut1)+"_cut2_"+str(bdt_cut2)+extra_name+" >> result_combine_3bin"+extra_name+".txt")


  # 4 Bin analysis 
  if run_4bin:
    os.system("rm result_combine_4bin"+extra_name+".txt")
    for icut1 in range(start_bin,8):
       for icut2 in range(icut1+1,9):
          for icut3 in range(icut2+1,10):
            # 4 Bin analysis
            datacard_4bin="imax    4 number of bins\njmax    * number of processes minus 1\nkmax    * number of nuisance parameters\n"
            datacard_4bin+="----------------------------\n---------------------------\n"
            datacard_4bin+="bin signal_region1 signal_region2 signal_region3 signal_region4\n"
            bdt_cut1 = icut1*step
            bdt_cut2 = icut2*step
            bdt_cut3 = icut3*step
            sgn1 = Getting_MC_exp(path_emu,nden_emu,sf_emu,start_val,bdt_cut1,BFxLxSigma) 
            bkg1= Getting_Bks_inSR(path_data_emu,"bin4a",start_val,icut1*0.1)
            sgn2 = Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut1,bdt_cut2,BFxLxSigma)    
            bkg2 = Getting_Bks_inSR(path_data_emu,"bin4b",bdt_cut1,bdt_cut2)
            sgn3 = Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut2,bdt_cut3,BFxLxSigma)
            bkg3 = Getting_Bks_inSR(path_data_emu,"bin4c",bdt_cut2,bdt_cut3)
            sgn4 = Getting_MC_exp(path_emu,nden_emu,sf_emu,bdt_cut3,1.0,BFxLxSigma)
            bkg4 = Getting_Bks_inSR(path_data_emu,"bin4d",bdt_cut3,1.0)
            datacard_4bin+="observation "+str(bkg1)+" "+str(bkg2)+" "+str(bkg3)+" "+str(bkg4)+"\n"
            datacard_4bin+="bin signal_region1 signal_region1 signal_region2 signal_region2 signal_region3 signal_region3 signal_region4 signal_region4\n"
            datacard_4bin+="process background signal background signal background signal background signal\n"
            datacard_4bin+="process 1 0 1 0 1 0 1 0\n"
            datacard_4bin+="rate "+str(bkg1)+" "+str(sgn1)+"  "+str(bkg2)+" "+str(sgn2)+" "+str(bkg3)+" "+str(sgn3)+" "+str(bkg4)+" "+str(sgn4)
            with open("combine_bin4_cut1_"+str(bdt_cut1)+"_cut2_"+str(bdt_cut2)+"_cut3_"+str(bdt_cut3)+extra_name+".txt",'w') as txt:
               txt.write(datacard_4bin)
               txt.close()
            os.system('echo bdt_cut '+str(bdt_cut1)+' bdt_cut2 '+str(bdt_cut2)+' bdt_cut3 '+str(bdt_cut3)+'>>result_combine_4bin'+extra_name+'.txt')
            os.system("./run_combine.sh combine_bin4_cut1_"+str(bdt_cut1)+"_cut2_"+str(bdt_cut2)+"_cut3_"+str(bdt_cut3)+extra_name+".txt combine_bin4_cut1_"+str(bdt_cut1)+"_cut2_"+str(bdt_cut2)+"_cut3_"+str(bdt_cut3)+extra_name+" >> result_combine_4bin"+extra_name+".txt")
  
