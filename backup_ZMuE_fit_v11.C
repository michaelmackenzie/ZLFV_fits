#include <fstream>
#include "fit_helper.h"

void DefaultVals(std::vector<RooRealVar> vars, std::vector<double> vals){
  for (int i =0; i< vars.size(); i++)
    vars[i].setVal(vals[i]);
}



std::pair<double,double> SignalSystematicsMaxMeanWidth(TString syst_file, TString cuts, TString syst_name, float min_fit_range, float max_fit_range, int nbin_data){
  
    TTree * syst_tree = get_tree("mytreefit",syst_file,cuts); 
    RooRealVar dilep_mass_syst("mass_ll_"+syst_name,"m(e,#mu)", (min_fit_range-max_fit_range)/2., min_fit_range , max_fit_range, "GeV/c^{2}");
    RooDataSet syst_dataset("syst_"+syst_name+"_dataset","",RooArgSet(dilep_mass_syst),RooFit::Import(*syst_tree));
    RooRealVar syst_cb_mean("sgn_cb_mean_"+syst_name,"",91.0e+00, 84.0e+00, 100.0e+00);
    RooRealVar syst_cb_width("syst_cb_width_"+syst_name,"",5., 0.1, 10.);
    RooRealVar syst_cb_a1("syst_cb_a1_"+syst_name,"",1.0, 0.1, 100.0);
    RooRealVar syst_cb_n1("syst_cb_n1_"+syst_name,"",1.0, 0.00001, 100.0);
    RooRealVar syst_cb_a2("syst_cb_a2_"+syst_name,"",1.0, 0.1, 100.0);
    RooRealVar syst_cb_n2("syst_cb_n2_"+syst_name,"",1.0, 0.1, 100.0);

    RooDoubleCB syst_sgn_PDF("sgn_PDF_"+syst_name+"_syst","cb",dilep_mass_syst,syst_cb_mean,syst_cb_width,syst_cb_a1,syst_cb_n1,syst_cb_a2,syst_cb_n2);
    RooRealVar n_syst("n_syst"+syst_name, "",1000,0,10000000000);
    RooAddPdf syst_esgn_PDF("esgn_PDF_"+syst_name,"",RooArgList(syst_sgn_PDF),RooArgList(n_syst));
    RooFitResult * syst_result = syst_esgn_PDF.fitTo(syst_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range , max_fit_range),RooFit::PrintLevel(-1));
    double mean=syst_cb_mean.getVal(), width=syst_cb_width.getVal();
    cout<<"signal systematic "+syst_name<<" mean "<<mean<<" width "<<width<<endl;
    auto syst_frame = dilep_mass_syst.frame();
    print_details (syst_result, syst_frame);
    syst_dataset.plotOn(syst_frame,RooFit::Binning(nbin_data));
    syst_esgn_PDF.plotOn(syst_frame,RooFit::LineColor(kBlue));
    save_plot(syst_frame,"m(#mu,e)","v11_prmtv_syst_"+syst_name);
    delete syst_result;
    return std::make_pair(mean,width);
}


RooDoubleCB SignalSystematicsFunction(TString varname, double max_shift_muon, double max_shift_electron,RooRealVar dilep_mass_out, RooRealVar sgn_cb_mean, RooRealVar sgn_cb_width, RooRealVar sgn_cb_a1, RooRealVar sgn_cb_n1, RooRealVar sgn_cb_a2, RooRealVar sgn_cb_n2){
   RooRealVar *mu_scale = new RooRealVar("nuisance_mu_scale","",0,-5,5);
   RooRealVar *ele_scale = new RooRealVar("nuisance_ele_scale","",0,-5,5);
   mu_scale->setConstant(true);
   ele_scale->setConstant(true);
   RooFormulaVar *mean_formula=new RooFormulaVar("mean_formula_"+varname, "", "@0*(1+"+TString(std::to_string(max_shift_muon))+"*@1)*(1+"+TString(std::to_string(max_shift_electron))+"*@2)", RooArgList(sgn_cb_mean,*mu_scale,*ele_scale));
   RooDoubleCB syst_sgn_PDF_out ("syst_sgn_PDF_"+varname,"cb",dilep_mass_out,*mean_formula,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
   return syst_sgn_PDF_out;
}


std::vector<RooRealVar> ChebParams(int order, TString varname){
   std::vector<RooRealVar> bkg_chebs;
   for (int i =0; i<order; i++){
     float def=0,min=-1.0,max=1.0;
     if (i==0)
        def=-1,min=-5.0,max=5.0;
     else if (i==1) def=0.6;
     else if (i==2) def=-0.2;
     RooRealVar bkg_cheb_x("bkg_cheb"+TString(std::to_string(order))+"_x"+TString(std::to_string(i))+"_"+varname, "", def,min,max);
     bkg_chebs.push_back(bkg_cheb_x);
   }
   return  bkg_chebs;
}

std::vector<RooRealVar> BstParams(int order, TString varname){
   std::vector<RooRealVar> bkg_bsts;
   for (int i =0; i<order; i++){
     RooRealVar bkg_bst_x("bkg_bst"+TString(std::to_string(order))+"_x"+TString(std::to_string(i))+"_"+varname, "", 1./TMath::Power(10,order), -25, 25);
     bkg_bsts.push_back(bkg_bst_x);
   }
   return  bkg_bsts;
}

void PlotFunctions(std::vector<RooAbsPdf*> pdfs, RooPlot * xframe, RooRealVar dilep_mass, std::vector<TString> names, std::vector<TString> legs, TLegend *leg, RooDataSet dataset, TString name, int nbin_data){
   for (int i=0; i<pdfs.size(); i++){
     pdfs[i]->plotOn(xframe,RooFit::LineColor(i+1),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name(names[i]));
     save_pull(xframe, dilep_mass, "m(e,#mu)","v11_"+names[i]);
     leg->AddEntry(xframe->findObject(names[i]),legs[i]);
   }
   dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::CutRange("left,right"),RooFit::Name("data"));
   leg->AddEntry(xframe->findObject("data"),"Data");
   save_plot(xframe,"m(e,#mu)","v11_sr_total_"+name,leg);
}

void FitFunctions(std::vector<RooAbsPdf*> pdfs, RooDataSet &dataset, RooPlot * xframe, vector<TString> names){
  for (int i=0; i<pdfs.size(); i++){
    RooFitResult * fit_result = pdfs[i]->fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));
    cout<<names[i]<<endl;
    print_details (fit_result, xframe); 
  }
}

std::pair<double,double> yield_calc( float nYld_total, RooRealVar dilep_mass, RooAbsPdf *pdf){
   RooAbsReal* rsgn = pdf->createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr"));
   double nYld_inSR = rsgn->getVal()*nYld_total;
   return make_pair(nYld_total,nYld_inSR);
}

int ZMuE_fit_v11(TString name="test", 
    TString data_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root",
    TString sgn_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root",
    TString sgn_mu_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEUP_bdt_v7_signal_mcRun1*.root", 
    TString sgn_mu_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEDOWN_bdt_v7_signal_mcRun1*.root", 
    TString sgn_ele_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEUP_bdt_v7_signal_mcRun1*.root", 
    TString sgn_ele_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEDOWN_bdt_v7_signal_mcRun1*.root", 
    TString xgbmin="0.9",TString xgbmax="1.0", bool plot_primitives=true, bool unblind=false, bool create_dc_input=false, float expected_Nsgn=350, TString outvar="mass_ll", bool syst_sgn=false, bool altfit_bkg=true, TString varname="bin", int ntoys=-1, float r_toys=1.0){

   gROOT->SetBatch(true);
   TString dilep_var_name="mass_ll";
   TString dilep_muup_name="mass_ll_Muon_scale_up";
   TString dilep_mudown_name="mass_ll_Muon_scale_down";
   TString dilep_eleup_name="mass_ll_Electron_scale_up";
   TString dilep_eledown_name="mass_ll_Electron_scale_down";
   std::vector<TString> syst_names{"Muon_scale_up", "Muon_scale_down", "Electron_scale_up", "Electron_scale_down"};
   std::vector<TString> syst_files{sgn_mu_up_file, sgn_mu_down_file, sgn_ele_up_file, sgn_ele_down_file};
  

   bool altfit_bkg_bst3=false;
   bool altfit_bkg_bst4=false;
   bool altfit_bkg_cheb4=false;
   bool altfit_bkg_cheb5=false;
   bool altfit_bkg_gamma=false;
   bool altfit_bkg_exp=true; 
 
   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron";
   TString tree_name="mytreefit";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;
   double blind_min = 86.;
   double blind_max = 96.;
   int nbin_data = 80;
   bool Verbose=false;
   bool Save_gen_plots=true;
   
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   
   // read trees
   TTree * sgn_tree = get_tree("mytreefit",sgn_file,cuts);
   TTree * data_tree = get_tree("mytreefit",data_file,cuts);
   TH1F* hdata = new TH1F("hdata","", nbin_data,min_fit_range , max_fit_range); 
   data_tree->Draw("mass_ll>>hdata");

  
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooDataSet sgn_dataset("sgn_dataset","sgn_dataset",RooArgSet(dilep_mass),RooFit::Import(*sgn_tree));
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");


   //////////////////////////// signal /////////////////////////////////////////
   RooWorkspace ws ("ws");
   ws.import(dilep_mass);
    
   RooRealVar sgn_cb_mean("sgn_cb_mean_"+varname,"",91.0e+00, 84.0e+00, 100.0e+00);
   RooRealVar sgn_cb_width("sgn_cb_width_"+varname,"",5., 0.1, 10.);
   RooRealVar sgn_cb_a1("sgn_cb_a1_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n1("sgn_cb_n1_"+varname,"",1.0, 0.00001, 100.0);
   RooRealVar sgn_cb_a2("sgn_cb_a2_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n2("sgn_cb_n2_"+varname,"",1.0, 0.1, 100.0);

   RooDoubleCB sgn_PDF("sgn_PDF","cb",dilep_mass,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2); 

   RooRealVar n_sgn("sgn_PDF_norm", "",1000,0,10e10);
   RooAddPdf esgn_PDF("sgn","esgnPDF",RooArgList(sgn_PDF),RooArgList(n_sgn));

   RooFitResult * sgn_result = esgn_PDF.fitTo(sgn_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range, max_fit_range),RooFit::PrintLevel(-1));
   sgn_cb_mean.setConstant(true);   sgn_cb_width.setConstant(true);
   sgn_cb_a1.setConstant(true);     sgn_cb_n1.setConstant(true);
   sgn_cb_a2.setConstant(true);     sgn_cb_n2.setConstant(true);
   RooDoubleCB sgn_PDF_out("sgn_PDF_"+varname,"cb",dilep_mass_out,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);

   auto sgn_frame = dilep_mass.frame();   
   print_details(sgn_result, sgn_frame);
   sgn_dataset.plotOn(sgn_frame,RooFit::Binning(nbin_data));
   esgn_PDF.plotOn(sgn_frame,RooFit::LineColor(kBlue));

   n_sgn.setVal(expected_Nsgn);
   if (plot_primitives)
     save_plot(sgn_frame,"m(#mu,e)","v11_prmtv_sgn_"+name);
   

   // signal systematics
   double mean_shift_mu_up=0,width_shift_mu_up=0,mean_shift_mu_down=0,width_shift_mu_down=0;
   double mean_shift_ele_up=0,width_shift_ele_up=0,mean_shift_ele_down=0,width_shift_ele_down=0;
   
   if (syst_sgn){
     for (int isyst=0; isyst<syst_files.size(); isyst++){
        std::pair<double,double> max_mean_width = SignalSystematicsMaxMeanWidth(syst_files[isyst], cuts, syst_names[isyst], min_fit_range, max_fit_range, nbin_data) ;
        if (isyst==0){ mean_shift_mu_up=max_mean_width.first; width_shift_mu_up=max_mean_width.second;}
        if (isyst==1){ mean_shift_mu_down=max_mean_width.first; width_shift_mu_down=max_mean_width.second;}
        if (isyst==2){ mean_shift_ele_up=max_mean_width.first; width_shift_ele_up=max_mean_width.second;}
        if (isyst==3){ mean_shift_ele_down=max_mean_width.first; width_shift_ele_down=max_mean_width.second;}
     }
   }

   
   if (create_dc_input){
      RooRealVar n_sgn_out("sgn_PDF_"+varname+"_norm", "",expected_Nsgn,0,10*expected_Nsgn);
      n_sgn_out.setConstant(true);
     
      RooWorkspace *wspace_sgn = new RooWorkspace("workspace_signal","workspace_signal");
      wspace_sgn->import(sgn_PDF_out);
      wspace_sgn->import(n_sgn_out);

     double max_shift_muon=0;
     double max_shift_electron=0;
     if (syst_sgn){
        max_shift_muon= fabs(sgn_cb_mean.getVal()-mean_shift_mu_up)>fabs(sgn_cb_mean.getVal()-mean_shift_mu_down) ? fabs(sgn_cb_mean.getVal()-mean_shift_mu_up)/sgn_cb_mean.getVal() : fabs(sgn_cb_mean.getVal()-mean_shift_mu_down)/sgn_cb_mean.getVal();
        max_shift_electron= fabs(sgn_cb_mean.getVal()-mean_shift_ele_up)>fabs(sgn_cb_mean.getVal()-mean_shift_ele_down) ? fabs(sgn_cb_mean.getVal()-mean_shift_ele_up)/sgn_cb_mean.getVal() : fabs(sgn_cb_mean.getVal()-mean_shift_ele_down)/sgn_cb_mean.getVal();
        cout<<"delta/M muon "<<max_shift_muon<<" electron "<<max_shift_electron<<endl;  
     }
//     RooDoubleCB syst_sgn_PDF_out= SignalSystematicsFunction(varname, max_shift_muon, max_shift_electron, dilep_mass_out,sgn_cb_mean, sgn_cb_width, sgn_cb_a1, sgn_cb_n1, sgn_cb_a2, sgn_cb_n2);
     RooRealVar *mu_scale = new RooRealVar("nuisance_mu_scale","",0,-5,5);
     RooRealVar *ele_scale = new RooRealVar("nuisance_ele_scale","",0,-5,5);
     mu_scale->setConstant(true);
     ele_scale->setConstant(true);
     RooFormulaVar *mean_formula=new RooFormulaVar("mean_formula_"+varname, "", "@0*(1+"+TString(std::to_string(max_shift_muon))+"*@1)*(1+"+TString(std::to_string(max_shift_electron))+"*@2)", RooArgList(sgn_cb_mean,*mu_scale,*ele_scale));
     RooDoubleCB syst_sgn_PDF_out ("syst_sgn_PDF_"+varname,"cb",dilep_mass_out,*mean_formula,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
     RooRealVar n_syst_sgn_out("syst_sgn_PDF_"+varname+"_norm", "",expected_Nsgn,0,10000000000);
     if (syst_sgn){
        wspace_sgn->import(syst_sgn_PDF_out);
        wspace_sgn->import(n_syst_sgn_out);
      }
      wspace_sgn->writeToFile("workspace_v11_sgn_"+name+".root");
   }
   
  
   ////////////////////////////////// Bkg ////////////////////////////////////
   RooDataSet dataset("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));
   ws.import(dataset); 

   dilep_mass.setRange("left",min_fit_range, blind_min); //70,88
   dilep_mass.setRange("right",blind_max, max_fit_range); //94,110
   dilep_mass.setRange("sr",blind_min,blind_max); //88,94
   dilep_mass.setRange("full",min_fit_range, max_fit_range);

   std::vector<RooAbsPdf*> bkg_functions;
   std::vector<TString> bkg_fnc_names;
   std::vector<TString> bkg_fnc_legs;
   std::vector<RooRealVar*> bkg_ampl;

   auto xframe = dilep_mass.frame();
   dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::MarkerColor(0),RooFit::LineColor(0));
   
   /////////////////////////////// cheb 3rd order (main) //////////////////////
   std::vector<RooRealVar> bkg_cheb3_params = ChebParams(3, varname);
   RooChebychev cheb3_bkgPDF("cheb3_bkgPDF","",dilep_mass,RooArgList(bkg_cheb3_params[0],bkg_cheb3_params[1],bkg_cheb3_params[2]));
   RooRealVar cheb3_n_bkg("cheb3_bkgPDF_"+varname+"_norm","",dataset.numEntries(),0,2*dataset.numEntries());
   RooAddPdf cheb3_ebkgPDF("cheb3_ebkgPDF","",RooArgList(cheb3_bkgPDF),RooArgList(cheb3_n_bkg));
   RooChebychev cheb3_bkgPDF_out("cheb3_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb3_params[0],bkg_cheb3_params[1],bkg_cheb3_params[2]));
   bkg_functions.push_back(&cheb3_ebkgPDF); 
   bkg_fnc_names.push_back("cheb3_"+name);
   bkg_fnc_legs.push_back("Chebychev 3");
   bkg_ampl.push_back(&cheb3_n_bkg);

   /////////////////////////////// cheb 4th order /////////////////////////////
   std::vector<RooRealVar> bkg_cheb4_params = ChebParams(4, varname);
   RooChebychev cheb4_bkgPDF("cheb4_bkgPDF","",dilep_mass,RooArgList(bkg_cheb4_params[0],bkg_cheb4_params[1],bkg_cheb4_params[2],bkg_cheb4_params[3]));
   RooRealVar cheb4_n_bkg("cheb4_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf cheb4_ebkgPDF("cheb4_ebkgPDF","",RooArgList(cheb4_bkgPDF),RooArgList(cheb4_n_bkg));
   RooChebychev cheb4_bkgPDF_out("cheb4_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb4_params[0],bkg_cheb4_params[1],bkg_cheb4_params[2],bkg_cheb4_params[3]));


   /////////////////////////////// cheb 5th order /////////////////////////////
   std::vector<RooRealVar> bkg_cheb5_params = ChebParams(5, varname);
   RooChebychev cheb5_bkgPDF("cheb5_bkgPDF","",dilep_mass,RooArgList(bkg_cheb5_params[0],bkg_cheb5_params[1],bkg_cheb5_params[2],bkg_cheb5_params[3],bkg_cheb5_params[4]));
   RooRealVar cheb5_n_bkg("cheb5_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf cheb5_ebkgPDF("cheb5_ebkgPDF","",RooArgList(cheb5_bkgPDF),RooArgList(cheb5_n_bkg));
   RooChebychev cheb5_bkgPDF_out("cheb5_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb5_params[0],bkg_cheb5_params[1],bkg_cheb5_params[2],bkg_cheb5_params[3],bkg_cheb5_params[4]));

   /////////////////////////////// 3rd order B. ///////////////////////////////
   std::vector<RooRealVar> bkg_bst3_params= BstParams(3, varname);
   RooBernsteinFast<3> bst3_bkgPDF("bst3_bkgPDF_"+varname, "", dilep_mass, RooArgList(bkg_bst3_params[0],bkg_bst3_params[1],bkg_bst3_params[2]) );
   bst3_bkgPDF.protectSubRange(true);
   RooRealVar bst3_n_bkg("bst3_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());   
   RooAddPdf bst3_ebkgPDF("bst3_ebkgPDF","",RooArgList(bst3_bkgPDF),RooArgList(bst3_n_bkg));
   RooBernsteinFast<3> bst3_bkgPDF_out("bst3_bkgPDF_"+varname, "", dilep_mass_out, RooArgList(bkg_bst3_params[0],bkg_bst3_params[1],bkg_bst3_params[2]) );

   /////////////////////////////// 4th order B. ////////////////////////////////
   std::vector<RooRealVar> bkg_bst4_params= BstParams(4, varname);
   RooBernsteinFast<4> bst4_bkgPDF("bst4_bkgPDF", "", dilep_mass, RooArgList(bkg_bst4_params[0],bkg_bst4_params[1],bkg_bst4_params[2],bkg_bst4_params[3]) );
   bst4_bkgPDF.protectSubRange(true);
   RooRealVar bst4_n_bkg("bst4_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());   
   RooAddPdf bst4_ebkgPDF("bst4_ebkgPDF","",RooArgList(bst4_bkgPDF),RooArgList(bst4_n_bkg));
   RooBernsteinFast<4> bst4_bkgPDF_out("bst4_bkgPDF_"+varname, "", dilep_mass_out, RooArgList(bkg_bst4_params[0],bkg_bst4_params[1],bkg_bst4_params[2],bkg_bst4_params[3]) );

   /////////////////////////////// gamma ////////////////////////////////////
   RooRealVar bkg_gamma_a("bkg_gamma_g_"+varname, "",0.7, 0.01, 20.); //5., 0.01, 10.); // 60., 50., 190.
   RooRealVar bkg_gamma_b("bkg_gamma_b_"+varname, "", 26, 0.01, 1500.0); //10, 0.01, 1500.0); // 0.05, 0.01, 0.1
   RooRealVar bkg_gamma_mu("bkg_gamma_mu_"+varname, "", -100, -1000, 70.); //1, -1000, 70.);
   RooGamma gamma_bkgPDF("gamma_bkgPDF","", dilep_mass, bkg_gamma_a,bkg_gamma_b, bkg_gamma_mu);
   RooRealVar gamma_n_bkg("gamma_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf gamma_ebkgPDF("gamma_ebkgPDF","",RooArgList(gamma_bkgPDF),RooArgList(gamma_n_bkg));
   RooGamma gamma_bkgPDF_out("gamma_bkgPDF_"+varname,"", dilep_mass_out, bkg_gamma_a,bkg_gamma_b, bkg_gamma_mu);

   /////////////////////////////// pol exp //////////////////////////////////// 
//   RooRealVar bkg_polexp_x0("bkg_polexp_x0_"+varname, "", -1.8357e-05,-10., 10.);
//   RooRealVar bkg_polexp_x1("bkg_polexp_x1_"+varname, "", -9.0861e-02, -10., 10.);
   RooRealVar bkg_polexp_x0("bkg_polexp_x0_"+varname, "", -0.0001,-10., 10.);
   RooRealVar bkg_polexp_x1("bkg_polexp_x1_"+varname, "", 0.1, -10., 10.);
   RooRealVar bkg_polexp_x2("bkg_polexp_x2_"+varname, "", -0.1, -10., 10.);
//   RooRealVar bkg_polexp_c0("bkg_polexp_c0_"+varname, "", 1.e3, 0., 1.e6);
   RooRealVar bkg_polexp_c0("bkg_polexp_c0_"+varname, "", 0.3, 0., 1.0);
   RooRealVar bkg_polexp_c1("bkg_polexp_c1_"+varname, "", 0.3, 0., 1.0);

   RooExponential polexp_bkgPDF1("polexp_bkgPDF1","", dilep_mass,bkg_polexp_x0 );
   RooExponential polexp_bkgPDF2("polexp_bkgPDF2","", dilep_mass,bkg_polexp_x1 );
   RooExponential polexp_bkgPDF3("polexp_bkgPDF3","", dilep_mass,bkg_polexp_x2 );
   RooAddPdf polexp_bkgPDF("polexp_bkgPDF","",RooArgList(polexp_bkgPDF1,polexp_bkgPDF2,polexp_bkgPDF3),RooArgList(bkg_polexp_c0,bkg_polexp_c1));
   
   RooRealVar polexp_n_bkg("polexp_bkgPDF_"+varname+"_norm","",dataset.numEntries()/2.,0,1.5*dataset.numEntries());
   RooAddPdf polexp_ebkgPDF("polexp_ebkgPDF","",RooArgList(polexp_bkgPDF),RooArgList(polexp_n_bkg));
   
    RooExponential polexp_bkgPDF1_out("polexp_bkgPDF1_"+varname,"", dilep_mass_out,bkg_polexp_x0 );
    RooExponential polexp_bkgPDF2_out("polexp_bkgPDF2_"+varname,"", dilep_mass_out,bkg_polexp_x1 );
    RooAddPdf polexp_bkgPDF_out("polexp_bkgPDF_"+varname,"",RooArgList(polexp_bkgPDF1_out,polexp_bkgPDF2_out),RooArgList(bkg_polexp_c0));



   /////////////////////////////// alternative fits bkg ////////////////////////
   if (altfit_bkg){
      if (altfit_bkg_cheb4){ 
         bkg_functions.push_back(&cheb4_ebkgPDF); 
         bkg_fnc_names.push_back("cheb4_"+name);
         bkg_fnc_legs.push_back("Chebychev 4");
         bkg_ampl.push_back(&cheb4_n_bkg);
      }
      if (altfit_bkg_cheb5){
         bkg_functions.push_back(&cheb5_ebkgPDF); 
         bkg_fnc_names.push_back("cheb5_"+name);
         bkg_fnc_legs.push_back("Chebychev 5");
         bkg_ampl.push_back(&cheb5_n_bkg);
      }
      if (altfit_bkg_bst3){
         bkg_functions.push_back(&bst3_ebkgPDF); 
         bkg_fnc_names.push_back("bst3_"+name);
         bkg_fnc_legs.push_back("Bernstein 3");
         bkg_ampl.push_back(&bst3_n_bkg);
      }
      if (altfit_bkg_bst4){
         bkg_functions.push_back(&bst4_ebkgPDF); 
         bkg_fnc_names.push_back("bst4_"+name);
         bkg_fnc_legs.push_back("Bernstein 4");
         bkg_ampl.push_back(&bst4_n_bkg);
       }
      if (altfit_bkg_exp){ 
         bkg_functions.push_back(&polexp_ebkgPDF); 
         bkg_fnc_names.push_back("expo_"+name);
         bkg_fnc_legs.push_back("#Sigma Expo");
         bkg_ampl.push_back(&polexp_n_bkg);
      }
      if (altfit_bkg_gamma){
         bkg_functions.push_back(&gamma_ebkgPDF); 
         bkg_fnc_names.push_back("gamma_"+name);
         bkg_fnc_legs.push_back("Gamma");
         bkg_ampl.push_back(&gamma_n_bkg);
      }
   }  
   //cout<<"before "<<bkg_cheb3_params[0].getVal()<<"  "<<bkg_cheb3_params[1].getVal()<<"  "<<bkg_cheb3_params[2].getVal()<<endl;
   FitFunctions(bkg_functions,dataset,xframe,bkg_fnc_legs);
   std::vector<double> final_params_cheb3 ={bkg_cheb3_params[0].getVal(),bkg_cheb3_params[1].getVal(),bkg_cheb3_params[2].getVal()};
   cout<<"after "<<final_params_cheb3[0]<<"  "<<final_params_cheb3[1]<<"  "<<final_params_cheb3[2]<<endl;
   //cout<<"after "<<bkg_cheb3_params[0].getVal()<<"  "<<bkg_cheb3_params[1].getVal()<<"  "<<bkg_cheb3_params[2].getVal()<<endl;

   /////////////////////////////// plotting /////////////////////////////////  
   TLegend * leg = new TLegend(0.6,0.6,0.9,0.85);
   PlotFunctions(bkg_functions, xframe,dilep_mass,bkg_fnc_names,bkg_fnc_legs,leg, dataset, name, nbin_data);
   esgn_PDF.plotOn(xframe,RooFit::LineColor(kBlue),RooFit::Normalization(1, RooAbsReal::RelativeExpected),RooFit::Name("sgn"));
   leg->AddEntry(xframe->findObject("sgn"),"Signal(DCB)");
   save_plot(xframe,"m(e,#mu)","v11_data_fit_and_sgn_"+name,leg);

   ///////////////////////////// output combine //////////////////////////////
   if (create_dc_input){
      RooDataSet * data_obs = cheb3_bkgPDF_out.generate(RooArgSet(dilep_mass_out),dataset.numEntries());
      RooDataSet * signal_obs_r1 = sgn_PDF_out.generate(RooArgSet(dilep_mass_out),expected_Nsgn);
      RooDataSet * signal_obs_r2 = sgn_PDF_out.generate(RooArgSet(dilep_mass_out),expected_Nsgn*2);
      RooDataSet * signal_obs_r5 = sgn_PDF_out.generate(RooArgSet(dilep_mass_out),expected_Nsgn*5);
      RooDataSet *data_obs_r1 =(RooDataSet*) data_obs->Clone("sim_data_obs_r1_"+varname);
      RooDataSet *data_obs_r2 =(RooDataSet*) data_obs->Clone("sim_data_obs_r2_"+varname);
      RooDataSet *data_obs_r5 =(RooDataSet*) data_obs->Clone("sim_data_obs_r5_"+varname);
      data_obs_r1->append( *signal_obs_r1);
      data_obs_r2->append( *signal_obs_r2);
      data_obs_r5->append( *signal_obs_r5);

      data_obs->SetName("sim_data_obs_"+varname);
      auto gen_frame = dilep_mass_out.frame();
      auto gen_frame2 = dilep_mass_out.frame();
      data_obs->plotOn(gen_frame,RooFit::Binning(40),RooFit::MarkerColor(1),RooFit::LineColor(1));
      save_plot(gen_frame,"m(e,#mu)","v11_genbkg_"+name);
      data_obs_r1->plotOn(gen_frame,RooFit::Binning(40),RooFit::MarkerColor(2),RooFit::LineColor(2));
      data_obs_r2->plotOn(gen_frame,RooFit::Binning(40),RooFit::MarkerColor(3),RooFit::LineColor(3));
      data_obs_r5->plotOn(gen_frame,RooFit::Binning(40),RooFit::MarkerColor(4),RooFit::LineColor(4));
      signal_obs_r1->plotOn(gen_frame2,RooFit::Binning(40),RooFit::MarkerColor(2),RooFit::LineColor(2));
      signal_obs_r2->plotOn(gen_frame2,RooFit::Binning(40),RooFit::MarkerColor(3),RooFit::LineColor(3));
      signal_obs_r5->plotOn(gen_frame2,RooFit::Binning(40),RooFit::MarkerColor(4),RooFit::LineColor(4));
      
      save_plot(gen_frame,"m(e,#mu)","v11_gentotal_"+name);
      save_plot(gen_frame2,"m(e,#mu)","v11_gensignal_"+name);
      dilep_mass_out.setBins((max_fit_range-min_fit_range)/0.25);
      RooDataHist hdata_sim("sim_binned_obs_"+varname,"",dilep_mass_out,*data_obs);
      hdata_sim.SetName("sim_binned_obs_"+varname);
     
      RooDataHist hdata_sim_r1("sim_binned_r1_obs_"+varname,"",dilep_mass_out,*data_obs_r1);
      hdata_sim_r1.SetName("sim_binned_r1_obs_"+varname); 
      RooArgList models_out;
      DefaultVals(bkg_cheb3_params,{-1,0.6,-0.2});
      models_out.add(cheb3_bkgPDF_out);

      if (altfit_bkg){
         if (altfit_bkg_cheb4){
            DefaultVals(bkg_cheb4_params,{-1,0.6,-0.2,0}); 
            models_out.add(cheb4_bkgPDF_out);
         }
         if (altfit_bkg_cheb5){
            DefaultVals(bkg_cheb5_params,{-1,0.6,-0.2,0,0});
            models_out.add(cheb5_bkgPDF_out);
         }
         if (altfit_bkg_bst3){
            DefaultVals(bkg_bst3_params,{0.001,0.001,0.001});
            models_out.add(bst3_bkgPDF_out);
         }
         if (altfit_bkg_bst4){
            DefaultVals(bkg_bst4_params,{0.0001,0.0001,0.0001,0.0001});  
            models_out.add(bst4_bkgPDF_out);
         }
         if (altfit_bkg_exp){
            DefaultVals({bkg_polexp_x0,bkg_polexp_x1,bkg_polexp_c0},{-1.8357e-05,-9.0861e-02,0.5});
            models_out.add(polexp_bkgPDF_out);
         }
         if (altfit_bkg_gamma){
            DefaultVals({bkg_gamma_a,bkg_gamma_b,bkg_gamma_mu},{0.7,26,-100});
            models_out.add(gamma_bkgPDF_out);
         }
      }
      RooCategory cat("pdfindex_"+varname, "");
      RooMultiPdf multipdf("multipdf_"+varname, "", cat, models_out);
      RooRealVar norm_out("multipdf_"+varname+"_norm","",dataset.numEntries(),0,2*dataset.numEntries());      
      RooWorkspace *wspace_bkg = new RooWorkspace("workspace_background","workspace_background");
       
      wspace_bkg->import(cat, RooFit::RecycleConflictNodes());
      wspace_bkg->import(multipdf, RooFit::RecycleConflictNodes());
      wspace_bkg->import(norm_out);
      wspace_bkg->import(cheb3_bkgPDF_out);
      wspace_bkg->import(cheb3_n_bkg);   
      wspace_bkg->import(polexp_bkgPDF_out);
      wspace_bkg->import(polexp_n_bkg);  

      wspace_bkg->import(*data_obs);
      wspace_bkg->import(*data_obs_r1);
      wspace_bkg->import(*data_obs_r2);
      wspace_bkg->import(*data_obs_r5);
      wspace_bkg->import(hdata_sim);
      wspace_bkg->import(hdata_sim_r1);
      wspace_bkg->writeToFile("workspace_v11_bkg_"+name+".root");

   }
   
 
   cout<<"\n\nYield for "+name<<endl<<"whole range    |    85-95 only"<<endl;
   cout<<"Total Evt count: "<<dataset.numEntries()<<endl;
   std::pair<double,double> nSgn = yield_calc( n_sgn.getVal(), dilep_mass, &sgn_PDF);
   cout<<" - nSgn "<<nSgn.first<<"   |  "<<nSgn.second<<endl;
   for (int i=0; i<bkg_functions.size(); i++){
      std::pair<double,double> nBkg = yield_calc( bkg_ampl[i]->getVal(), dilep_mass, bkg_functions[i]);
      cout<<" - nBkg("+bkg_fnc_legs[i]+") "<<nBkg.first<<"  |  "<<nBkg.second<<endl;
   }
   cout<<" -- Chebychev3 Bkg main "<<endl;
   


   if (ntoys>0){
     cout<<"pull test"<<endl;
     TH1F * htoy_cheb3_nevt = new TH1F("htoy_cheb3_nevt","",100,dataset.numEntries()-5.*TMath::Sqrt(dataset.numEntries()),dataset.numEntries()+5.*TMath::Sqrt(dataset.numEntries()));
     RooRealVar toys_n_gen("toys_n_gen","",expected_Nsgn*r_toys);
     RooAddPdf total_pdf_toy("total_pdf_toy","",RooArgList(polexp_bkgPDF,sgn_PDF),RooArgList(polexp_n_bkg,toys_n_gen));
     cout<<toys_n_gen.getVal()*r_toys<<endl;
//     DefaultVals(bkg_cheb4_params,final_params_cheb3);
    //      DefaultVals(bkg_cheb3_params,{1,5,10});
  //   std::vector<RooRealVar> bkg_cheb3_params2 = ChebParams(3, varname+"2");
//     RooChebychev cheb3_bkgPDF_out2("cheb3_bkgPDF2_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb3_params2[0],bkg_cheb3_params2[1],bkg_cheb3_params2[2]));
     
     auto toy_cheb3_datasets = generate_pseudo_data(ntoys, dataset.numEntries(), total_pdf_toy, dilep_mass, true );
          

     for (int itoy=0; itoy<toy_cheb3_datasets.size(); itoy++){
       htoy_cheb3_nevt->Fill(toy_cheb3_datasets[itoy].numEntries());    
       if (Save_gen_plots){
          auto gendebug_frame = dilep_mass.frame();
          toy_cheb3_datasets[itoy].plotOn(gendebug_frame);
          TCanvas * cgentoys_debug = new TCanvas("cgentoys_debug_"+TString(std::to_string(itoy)),"",700,700);
          gendebug_frame->Draw();
          cgentoys_debug->SaveAs("cgentoys_debug_"+TString(std::to_string(itoy))+".png");
       }
     }

     TCanvas * ctoy_nevt = new TCanvas("ctoy_nevt","",700,700);
     htoy_cheb3_nevt->Draw();
     ctoy_nevt->SaveAs("ctoy_nevt.png");
    
     TH1F * hdiff_sgn0_exp = new TH1F("hdiff_sgn0_exp","",100,dataset.numEntries()*0.8,dataset.numEntries()*1.2);
     TH1F * hpull_sgn0_exp = new TH1F("hpull_sgn0_exp","",100,-5,5);
//     RooAddPdf total_pdf_exp_bkg("total_pdf_exp_bkg","",RooArgList(polexp_bkgPDF,sgn_PDF),RooArgList(polexp_n_bkg,n_sgn));
     
     int itoy=0;
     for (auto toy_cheb3: toy_cheb3_datasets ){
        if (altfit_bkg && altfit_bkg_exp){
     //     DefaultVals({bkg_polexp_x0,bkg_polexp_x1,bkg_polexp_c0,bkg_polexp_c1},{-1.8357e-05,-9.0861e-02,1.e3,1.e3});
//          polexp_ebkgPDF.fitTo(toy_cheb3,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"));
  //        total_pdf_exp_bkg.fitTo(toy_cheb3,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"));
         RooRealVar bkg_polexp2_x0("bkg_polexp2_x0_"+varname, "", -1.8357e-05,-10., 10.);
         RooRealVar bkg_polexp2_x1("bkg_polexp2_x1_"+varname, "", -9.0861e-02, -10., 10.);
         RooRealVar bkg_polexp2_c0("bkg_polexp2_c0_"+varname, "", 0.5, 0., 1.0);

         RooExponential polexp2_bkgPDF1("polexp2_bkgPDF1","", dilep_mass,bkg_polexp2_x0 );
         RooExponential polexp2_bkgPDF2("polexp2_bkgPDF2","", dilep_mass,bkg_polexp2_x1 );
         RooAddPdf polexp2_bkgPDF("polexp2_bkgPDF","",RooArgList(polexp2_bkgPDF1,polexp2_bkgPDF2),RooArgList(bkg_polexp2_c0));
         n_sgn.setVal(0);
         n_sgn.setRange(0,toys_n_gen.getVal()*2);         
         RooAddPdf total_pdf_exp_bkg("total_pdf_exp_bkg","",RooArgList(cheb3_bkgPDF,sgn_PDF),RooArgList(cheb3_n_bkg,n_sgn));
         auto fitdebug_frame = dilep_mass.frame();
         total_pdf_exp_bkg.fitTo(toy_cheb3,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"));
         toy_cheb3.plotOn(fitdebug_frame);
         total_pdf_exp_bkg.plotOn(fitdebug_frame);
         TCanvas * cfittoys_debug = new TCanvas("cfittoys_debug_"+TString(std::to_string(itoy)),"",700,700);
          fitdebug_frame->Draw();
          cfittoys_debug->SaveAs("cfittoys_debug_"+TString(std::to_string(itoy))+".png");
//          hdiff_sgn0_exp->Fill(polexp_n_bkg.getVal()-toy_cheb3.numEntries());

          hpull_sgn0_exp->Fill((n_sgn.getVal()-toys_n_gen.getVal())/n_sgn.getError());
          cout<<n_sgn.getVal()<<" - "<<toys_n_gen.getVal()<<" / "<<n_sgn.getError()<<" = "<<(n_sgn.getVal()-toys_n_gen.getVal())/n_sgn.getError()<<endl;
          itoy+=1;
        }  
     }
      TCanvas * cdiff_bkgonly_exp_cheb3 = new TCanvas("cdiff_bkgonly_exp_cheb3","",700,700);
      hpull_sgn0_exp->Draw();
      cdiff_bkgonly_exp_cheb3->SaveAs("cpull.png");
      
/*     TCanvas * cdiff_bkgonly_exp_cheb3 = new TCanvas("cdiff_bkgonly_exp_cheb3","",700,700);
     hdiff_bkgonly_exp->Draw();
     cdiff_bkgonly_exp_cheb3->SaveAs("cdiff_bkgonly_exp_cheb3.png");
     TCanvas * cpull_sgn0_exp_cheb3 = new TCanvas("cpull_sgn0_exp_cheb3","",700,700);
     hpull_sgn0_exp->Draw();
     cpull_sgn0_exp_cheb3->SaveAs("cpull_sgn0_exp_cheb3.png");*/

   }
 return 0;
} 
