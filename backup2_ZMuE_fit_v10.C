#include <fstream>
#include "fit_helper.h"






int ZMuE_fit_v10(TString name="test",  TString data_file="",TString sgn_file="",TString sgn_mu_up_file="", TString sgn_mu_down_file="", TString sgn_ele_up_file="", TString sgn_ele_down_file="", TString xgbmin="0.75",TString xgbmax="1.01", bool plot_primitives=true, bool unblind=false, bool create_dc_input=false, float expected_Nsgn=350, TString outvar="mass_ll", bool syst_sgn=true, bool altfit_bkg=true, TString varname="bin", int ntoys=0){

   gROOT->SetBatch(true);
   TString dilep_var_name="mass_ll";
   TString dilep_muup_name="mass_ll_Muon_scale_up";
   TString dilep_mudown_name="mass_ll_Muon_scale_down";
   TString dilep_eleup_name="mass_ll_Electron_scale_up";
   TString dilep_eledown_name="mass_ll_Electron_scale_down";
   bool altfit_bkg_bst3=true;
   bool altfit_bkg_bst4=false;
   bool altfit_bkg_cheb4=false;
   bool altfit_bkg_cheb5=false;
   bool altfit_bkg_gamma=false;
   bool altfit_bkg_pol4=false; 
   bool altfit_bkg_exp=false;
   
 
   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron";
   TString tree_name="mytreefit";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;
   double blind_min = 86.;
   double blind_max = 96.;
   int nbin_data = 80;
   bool Verbose=false;
   
   

////////////////////////////////////////////////////////////////////////////////

   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   
   // read trees
   TTree * sgn_tree = get_tree("mytreefit",sgn_file,cuts);
   TTree * data_tree = get_tree("mytreefit",data_file,cuts);
   TH1F* hdata = new TH1F("hdata","",80,70,110); 
   data_tree->Draw("mass_ll>>hdata");

  
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooDataSet sgn_dataset("sgn_dataset","sgn_dataset",RooArgSet(dilep_mass),RooFit::Import(*sgn_tree));
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
//   RooDataHist* datasetHist = new RooDataHist("blind_data_hist",  "Blind Data Hist", RooArgList(dilep_mass), hdata);

   //////////////////////////// signal ////////////////////////////// 
   RooWorkspace ws ("ws");
   ws.import(dilep_mass);
    
   RooRealVar sgn_cb_mean("sgn_cb_mean_"+varname,"",91.0e+00, 84.0e+00, 100.0e+00);
   RooRealVar sgn_cb_width("sgn_cb_width_"+varname,"",5., 0.1, 10.);
   RooRealVar sgn_cb_a1("sgn_cb_a1_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n1("sgn_cb_n1_"+varname,"",1.0, 0.00001, 100.0);
   RooRealVar sgn_cb_a2("sgn_cb_a2_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n2("sgn_cb_n2_"+varname,"",1.0, 0.1, 100.0);

   RooDoubleCB sgn_PDF("sgn_PDF","cb",dilep_mass,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2); 

   RooRealVar n_sgn("n_sgn", "n_sgn",1000,0,10000000000);
   RooAddPdf esgn_PDF("sgn","esgnPDF",RooArgList(sgn_PDF),RooArgList(n_sgn));


   // signal systematics
   RooFitResult * syst_result;
   std::vector<TString> syst_names{"Muon_scale_up", "Muon_scale_down", "Electron_scale_up", "Electron_scale_down"};
   std::vector<TString> syst_files{sgn_mu_up_file, sgn_mu_down_file, sgn_ele_up_file, sgn_ele_down_file};

   double mean_shift_mu_up=0,width_shift_mu_up=0,mean_shift_mu_down=0,width_shift_mu_down=0;
   double mean_shift_ele_up=0,width_shift_ele_up=0,mean_shift_ele_down=0,width_shift_ele_down=0;
   
   if (syst_sgn){
     for (int isyst=0; isyst<syst_files.size(); isyst++){
        TTree * syst_tree = get_tree("mytreefit",syst_files[isyst],cuts); 
        RooRealVar dilep_mass_syst("mass_ll_"+syst_names[isyst],"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
        RooDataSet syst_dataset("syst_"+syst_names[isyst]+"_dataset","",RooArgSet(dilep_mass_syst),RooFit::Import(*syst_tree));
        RooDoubleCB syst_sgn_PDF("sgn_PDF_"+syst_names[isyst]+"_syst","cb",dilep_mass_syst,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
        RooAddPdf syst_esgn_PDF("esgn_PDF_"+syst_names[isyst],"",RooArgList(syst_sgn_PDF),RooArgList(n_sgn));
        syst_result = syst_esgn_PDF.fitTo(syst_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(70,110),RooFit::PrintLevel(-1));
        cout<<"signal systematic "+syst_names[isyst]<<endl;
        auto syst_frame = dilep_mass_syst.frame();
        print_details (syst_result, syst_frame);
        syst_dataset.plotOn(syst_frame,RooFit::Binning(nbin_data));
        syst_esgn_PDF.plotOn(syst_frame,RooFit::LineColor(kBlue));
        save_plot(syst_frame,"m(#mu,e)","v10_prmtv_syst_"+syst_names[isyst]);
        cout<<" mean "<<sgn_cb_mean.getVal()<<"  width "<<sgn_cb_width.getVal()<<endl;            
        if (isyst==0){ mean_shift_mu_up=sgn_cb_mean.getVal(); width_shift_mu_up=sgn_cb_width.getVal();}
        if (isyst==1){ mean_shift_mu_down=sgn_cb_mean.getVal(); width_shift_mu_down=sgn_cb_width.getVal();}
        if (isyst==2){ mean_shift_ele_up=sgn_cb_mean.getVal(); width_shift_ele_up=sgn_cb_width.getVal();}
        if (isyst==3){ mean_shift_ele_down=sgn_cb_mean.getVal(); width_shift_ele_down=sgn_cb_width.getVal();}
     }
   }
   // nominal fit
   sgn_cb_mean.setVal(91.0);
   sgn_cb_width.setVal(5.);
   sgn_cb_a1.setVal(1.0);
   sgn_cb_n1.setVal(1.0);
   sgn_cb_a2.setVal(1.0);
   sgn_cb_n2.setVal(1.0);

   RooFitResult * sgn_result = esgn_PDF.fitTo(sgn_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(70,110),RooFit::PrintLevel(-1));

   auto sgn_frame = dilep_mass.frame();   
   print_details(sgn_result, sgn_frame);
   sgn_dataset.plotOn(sgn_frame,RooFit::Binning(nbin_data));
   esgn_PDF.plotOn(sgn_frame,RooFit::LineColor(kBlue));
   
   if (plot_primitives)
     save_plot(sgn_frame,"m(#mu,e)","v10_prmtv_sgn_"+name);
   n_sgn.setVal(expected_Nsgn);

   
   if (create_dc_input){
      sgn_cb_mean.setConstant(true);
      sgn_cb_width.setConstant(true);
      sgn_cb_a1.setConstant(true);
      sgn_cb_n1.setConstant(true);
      sgn_cb_a2.setConstant(true);
      sgn_cb_n2.setConstant(true);
      RooDoubleCB sgn_PDF_out("sgn_PDF_"+varname,"cb",dilep_mass_out,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
      RooRealVar n_sgn_out("sgn_PDF_"+varname+"_norm", "",expected_Nsgn,0,10000000000);
      n_sgn_out.setConstant(true);


      double max_shift_muon= fabs(sgn_cb_mean.getVal()-mean_shift_mu_up)>fabs(sgn_cb_mean.getVal()-mean_shift_mu_down) ? fabs(sgn_cb_mean.getVal()-mean_shift_mu_up)/sgn_cb_mean.getVal() : fabs(sgn_cb_mean.getVal()-mean_shift_mu_down)/sgn_cb_mean.getVal();
      double max_shift_electron= fabs(sgn_cb_mean.getVal()-mean_shift_ele_up)>fabs(sgn_cb_mean.getVal()-mean_shift_ele_down) ? fabs(sgn_cb_mean.getVal()-mean_shift_ele_up)/sgn_cb_mean.getVal() : fabs(sgn_cb_mean.getVal()-mean_shift_ele_down)/sgn_cb_mean.getVal();
      cout<<"delta/M muon "<<max_shift_muon<<" electron "<<max_shift_electron<<endl;  
      RooRealVar mu_scale("nuisance_mu_scale","",0,-5,5);
      RooRealVar ele_scale("nuisance_ele_scale","",0,-5,5);
      mu_scale.setConstant(true);
      ele_scale.setConstant(true);
      RooFormulaVar mean_formula("mean_formula_"+varname, "", "@0*(1+"+TString(std::to_string(max_shift_muon))+"*@1)*(1+"+TString(std::to_string(max_shift_electron))+"*@2)", RooArgList(sgn_cb_mean,mu_scale,ele_scale));

      RooDoubleCB syst_sgn_PDF_out("syst_sgn_PDF_"+varname,"cb",dilep_mass_out,mean_formula,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
      RooRealVar n_syst_sgn_out("syst_sgn_PDF_"+varname+"_norm", "",expected_Nsgn,0,10000000000);

      n_syst_sgn_out.setConstant(true);
      RooWorkspace *wspace_sgn = new RooWorkspace("workspace_signal","workspace_signal");
      wspace_sgn->import(sgn_PDF_out);
      wspace_sgn->import(n_sgn_out); 
      if (syst_sgn){
        wspace_sgn->import(syst_sgn_PDF_out);
        wspace_sgn->import(n_syst_sgn_out);
      }
      wspace_sgn->writeToFile("workspace_v10_sgn_"+name+".root");
   }
   
  
   ////////////////////////////////// Bkg ////////////////////////////////////
   RooDataSet dataset("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));
   ws.import(dataset); 

   dilep_mass.setRange("left",70, blind_min); //70,88
   dilep_mass.setRange("right",blind_max,110); //94,110
   dilep_mass.setRange("sr",blind_min,blind_max); //88,94
   dilep_mass.setRange("full",70,110);

   /////////////////////////////// cheb 3rd order (main) //////////////////////
   RooRealVar bkg_cheb3_x0("bkg_cheb3_x0_"+varname, "", -1 , -2.0, 2.0);
   RooRealVar bkg_cheb3_x1("bkg_cheb3_x1_"+varname, "", 0.6,  -1., 1.);
   RooRealVar bkg_cheb3_x2("bkg_cheb3_x2_"+varname, "", -0.2, -1., 1.);


   RooChebychev cheb3_bkgPDF("cheb3_bkgPDF","",dilep_mass,RooArgList(bkg_cheb3_x0,bkg_cheb3_x1,bkg_cheb3_x2));
   RooRealVar cheb3_n_bkg("cheb3_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());
   RooAddPdf cheb3_ebkgPDF("cheb3_ebkgPDF","",RooArgList(cheb3_bkgPDF),RooArgList(cheb3_n_bkg));
 
   RooFitResult * cheb3_result = cheb3_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));

   /////////////////////////////// cheb 4th order //////////////////////
   RooRealVar bkg_cheb4_x0("bkg_cheb4_x0_"+varname, "", -1 , -2.0, 2.0);
   RooRealVar bkg_cheb4_x1("bkg_cheb4_x1_"+varname, "", 0.6,  -1., 1.);
   RooRealVar bkg_cheb4_x2("bkg_cheb4_x2_"+varname, "", -0.2, -1., 1.);
   RooRealVar bkg_cheb4_x3("bkg_cheb4_x3_"+varname, "", 0, -1., 1.);
   RooRealVar bkg_cheb4_x4("bkg_cheb4_x4_"+varname, "", 0, -1., 1.);
   RooRealVar bkg_cheb4_x5("bkg_cheb4_x5_"+varname, "", 0, -1., 1.);

   RooChebychev cheb4_bkgPDF("cheb4_bkgPDF","",dilep_mass,RooArgList(bkg_cheb4_x0,bkg_cheb4_x1,bkg_cheb4_x2,bkg_cheb4_x3));
   RooRealVar cheb4_n_bkg("cheb4_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf cheb4_ebkgPDF("cheb4_ebkgPDF","",RooArgList(cheb4_bkgPDF),RooArgList(cheb4_n_bkg));
   RooFitResult * cheb4_result;
   if (altfit_bkg && altfit_bkg_cheb4)
      cheb4_result= cheb4_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));



   /////////////////////////////// cheb 5th order //////////////////////
   RooRealVar bkg_cheb5_x0("bkg_cheb5_x0_"+varname, "", -1 , -2.0, 2.0);
   RooRealVar bkg_cheb5_x1("bkg_cheb5_x1_"+varname, "", 0.6,  -1., 1.);
   RooRealVar bkg_cheb5_x2("bkg_cheb5_x2_"+varname, "", -0.2, -1., 1.);
   RooRealVar bkg_cheb5_x3("bkg_cheb5_x3_"+varname, "", 0, -1., 1.);
   RooRealVar bkg_cheb5_x4("bkg_cheb5_x4_"+varname, "", 0, -1., 1.);

   RooChebychev cheb5_bkgPDF("cheb5_bkgPDF","",dilep_mass,RooArgList(bkg_cheb5_x0,bkg_cheb5_x1,bkg_cheb5_x2,bkg_cheb5_x3,bkg_cheb5_x4));
   RooRealVar cheb5_n_bkg("cheb5_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf cheb5_ebkgPDF("cheb5_ebkgPDF","",RooArgList(cheb5_bkgPDF),RooArgList(cheb5_n_bkg));

   RooFitResult * cheb5_result;
   if (altfit_bkg && altfit_bkg_cheb5)
     cheb5_result = cheb5_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));   

   /////////////////////////////// 3rd order B. ////////////////////////////////
   RooRealVar bkg_bst3_x0("bkg_bst3_x0_"+varname, "bkg_bst3_x0",1./TMath::Power(10,3), -25, 25);
   RooRealVar bkg_bst3_x1("bkg_bst3_x1_"+varname, "bkg_bst3_x1",1./TMath::Power(10,4), -25, 25);
   RooRealVar bkg_bst3_x2("bkg_bst3_x2_"+varname, "bkg_bst3_x2",1./TMath::Power(10,4), -25, 25);
   RooBernsteinFast<3> bst3_bkgPDF("bst3_bkgPDF", "", dilep_mass, RooArgList(bkg_bst3_x0,bkg_bst3_x1,bkg_bst3_x2) );
   bst3_bkgPDF.protectSubRange(true);

   RooRealVar bst3_n_bkg("bst3_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());   
   RooAddPdf bst3_ebkgPDF("bst3_ebkgPDF","",RooArgList(bst3_bkgPDF),RooArgList(bst3_n_bkg));

   RooFitResult * bst3_result;
   if (altfit_bkg && altfit_bkg_bst3)
      bst3_result = bst3_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));



   /////////////////////////////// 4th order B. ////////////////////////////////
   RooRealVar bkg_bst4_x0("bkg_bst4_x0_"+varname, "bkg_bst4_x0",1./TMath::Power(10,4), -25, 25);
   RooRealVar bkg_bst4_x1("bkg_bst4_x1_"+varname, "bkg_bst4_x1",1./TMath::Power(10,4), -25, 25);
   RooRealVar bkg_bst4_x2("bkg_bst4_x2_"+varname, "bkg_bst4_x2",1./TMath::Power(10,4), -25, 25);
   RooRealVar bkg_bst4_x3("bkg_bst4_x3_"+varname, "bkg_bst4_x3",1./TMath::Power(10,4), -25, 25);
   RooBernsteinFast<4> bst4_bkgPDF("bst4_bkgPDF", "", dilep_mass, RooArgList(bkg_bst4_x0,bkg_bst4_x1,bkg_bst4_x2,bkg_bst4_x3) );
   bst4_bkgPDF.protectSubRange(true);

   RooRealVar bst4_n_bkg("bst4_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());   
   RooAddPdf bst4_ebkgPDF("bst4_ebkgPDF","",RooArgList(bst4_bkgPDF),RooArgList(bst4_n_bkg));

   RooFitResult * bst4_result;
   if (altfit_bkg && altfit_bkg_bst4)
      bst4_result = bst4_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));

   ////////////////////////// Polynomial ///////////////////////////
   RooRealVar bkg_pol4_x0("bkg_pol4_x0_"+varname, "", 1700000 ,-500., 100000000.);
   RooRealVar bkg_pol4_x1("bkg_pol4_x1_"+varname, "", -2000, -1000000., 10000000.);
   RooRealVar bkg_pol4_x2("bkg_pol4_x2_"+varname, "", 50.5, -1000., 1000.);
   RooRealVar bkg_pol4_x3("bkg_pol4_x3_"+varname, "", -0.5, -10., 1.);
   RooRealVar bkg_pol4_x4("bkg_pol4_x4_"+varname, "", 5.549e-1, -100., 0.);
   RooRealVar bkg_pol4_x5("bkg_pol4_x5_"+varname, "", 5.549e-1, -100., 100.);


   RooPolynomial pol4_bkgPDF("pol4_bkgPDF","",dilep_mass,RooArgList(bkg_pol4_x0,bkg_pol4_x1, bkg_pol4_x2,bkg_pol4_x3 ));
   RooRealVar pol4_n_bkg("pol4_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf pol4_ebkgPDF("pol4_ebkgPDF","",RooArgList(pol4_bkgPDF),RooArgList(pol4_n_bkg));

   RooFitResult * pol4_result;
   if (altfit_bkg && altfit_bkg_pol4)
      pol4_result = pol4_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));



   /////////////////////////////// gamma ////////////////////////////////////
   RooRealVar bkg_gamma_a("bkg_gamma_g_"+varname, "",0.7, 0.01, 20.); //5., 0.01, 10.); // 60., 50., 190.
   RooRealVar bkg_gamma_b("bkg_gamma_b_"+varname, "", 26, 0.01, 1500.0); //10, 0.01, 1500.0); // 0.05, 0.01, 0.1
   RooRealVar bkg_gamma_mu("bkg_gamma_mu_"+varname, "", -100, -1000, 70.); //1, -1000, 70.);

   RooGamma gamma_bkgPDF("gamma_bkgPDF","", dilep_mass, bkg_gamma_a,bkg_gamma_b, bkg_gamma_mu);
   RooRealVar gamma_n_bkg("gamma_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooAddPdf gamma_ebkgPDF("gamma_ebkgPDF","",RooArgList(gamma_bkgPDF),RooArgList(gamma_n_bkg));

   RooFitResult * gamma_result;
   if (altfit_bkg && altfit_bkg_gamma) 
      gamma_result = gamma_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"));

   /////////////////////////////// pol exp //////////////////////////////////// 
   RooRealVar bkg_polexp_x0("bkg_polexp_x0_"+varname, "", -1.8357e-05,-10., 10.);
   RooRealVar bkg_polexp_x1("bkg_polexp_x1_"+varname, "", -9.0861e-02, -10., 10.);
   RooRealVar bkg_polexp_c0("bkg_polexp_c0_"+varname, "", 1.e3, 0., 1.e6);
   RooRealVar bkg_polexp_c1("bkg_polexp_c1_"+varname, "", 1.e3, 0., 1.e6);

   RooExponential polexp_bkgPDF1("polexp_bkgPDF1","", dilep_mass,bkg_polexp_x0 );
   RooExponential polexp_bkgPDF2("polexp_bkgPDF2","", dilep_mass,bkg_polexp_x1 );
   RooAddPdf polexp_bkgPDF("polexp_bkgPDF","",RooArgList(polexp_bkgPDF1,polexp_bkgPDF2),RooArgList(bkg_polexp_c0,bkg_polexp_c1),false);
   
   RooRealVar polexp_n_bkg("polexp_bkgPDF_norm","",dataset.numEntries()/2.,0,1.5*dataset.numEntries());
   RooAddPdf polexp_ebkgPDF("polexp_ebkgPDF","",RooArgList(polexp_bkgPDF),RooArgList(polexp_n_bkg));
   polexp_bkgPDF.Print();
   RooFitResult * polexp_result;
   if (altfit_bkg && altfit_bkg_exp)
      polexp_result = polexp_ebkgPDF.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));


   /////////////////////////////// plotting /////////////////////////////////  
   auto xframe = dilep_mass.frame();
   TLegend * leg = new TLegend(0.6,0.6,0.9,0.85);

   dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::MarkerColor(0),RooFit::LineColor(0));
   cheb3_ebkgPDF.plotOn(xframe,RooFit::LineColor(kRed-5),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("cheb3"));
   save_pull(xframe, dilep_mass, "m(e,#mu)","v10_cheb3_"+name);
   leg->AddEntry(xframe->findObject("cheb3"),"Chebychev 3");
   print_details (cheb3_result, xframe);

   if (altfit_bkg){
     if (altfit_bkg_cheb4){
        cheb4_ebkgPDF.plotOn(xframe,RooFit::LineColor(kRed),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("cheb4"));
        save_pull(xframe, dilep_mass, "m(e,#mu)","v10_cheb4_"+name);
        leg->AddEntry(xframe->findObject("cheb4"),"Chebychev 4");
        print_details (cheb4_result, xframe);
     }
     if (altfit_bkg_cheb5){
        cheb5_ebkgPDF.plotOn(xframe,RooFit::LineColor(kRed+5),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("cheb5"));
        save_pull(xframe, dilep_mass, "m(e,#mu)","v10_cheb5_"+name);
        leg->AddEntry(xframe->findObject("cheb5"),"Chebychev 5");
        print_details (cheb5_result, xframe);
     }
     if (altfit_bkg_bst4){
        bst4_ebkgPDF.plotOn(xframe,RooFit::LineColor(kMagenta),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("bst4"));
        save_pull(xframe, dilep_mass, "m(e,#mu)","v10_bst4_"+name);
        print_details (bst4_result, xframe);
        leg->AddEntry(xframe->findObject("bst4"),"Bernstein 4");
     }
    if (altfit_bkg_bst3){
        bst3_ebkgPDF.plotOn(xframe,RooFit::LineColor(kMagenta-5),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("bst3"));
        save_pull(xframe, dilep_mass, "m(e,#mu)","v10_bst3_"+name);
        print_details (bst3_result, xframe);
        leg->AddEntry(xframe->findObject("bst3"),"Bernstein 3");
     }
     if (altfit_bkg_pol4){
         pol4_ebkgPDF.plotOn(xframe,RooFit::LineColor(kBlack),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("pol4"));
         save_pull(xframe, dilep_mass, "m(e,#mu)","v10_pol4_"+name);
         print_details (pol4_result, xframe);
          leg->AddEntry(xframe->findObject("pol4"),"Polynomial 4");
     }
     if (altfit_bkg_gamma){
        gamma_ebkgPDF.plotOn(xframe,RooFit::LineColor(kGreen+10),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("gamma"));
        save_pull(xframe, dilep_mass, "m(e,#mu)","v10_gamma_"+name);
        print_details (gamma_result, xframe);     
        leg->AddEntry(xframe->findObject("gamma"),"Gamma");
     }
     if (altfit_bkg_exp){
        polexp_ebkgPDF.plotOn(xframe,RooFit::LineColor(kGreen),RooFit::Range("full"),RooFit::NormRange("left,right"),RooFit::Name("polexp"));
        save_pull(xframe, dilep_mass, "m(e,#mu)","v10_polexp_"+name);
        print_details (polexp_result, xframe);
        leg->AddEntry(xframe->findObject("polexp"),"Exponential");
     }
   }
   
   dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::CutRange("left,right"),RooFit::Name("data"));
   leg->AddEntry(xframe->findObject("data"),"Data");

   save_plot(xframe,"m(e,#mu)","v10_sr_total_"+name,leg);

   if (create_dc_input){
      RooBernsteinFast<3> bst3_bkgPDF_out("bst3_bkgPDF_"+varname, "", dilep_mass_out, RooArgList(bkg_bst3_x0,bkg_bst3_x1,bkg_bst3_x2) );
      RooBernsteinFast<4> bst4_bkgPDF_out("bst4_bkgPDF_"+varname, "", dilep_mass_out, RooArgList(bkg_bst4_x0,bkg_bst4_x1,bkg_bst4_x2,bkg_bst4_x3) );
      RooChebychev cheb3_bkgPDF_out("cheb3_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb3_x0,bkg_cheb3_x1,bkg_cheb3_x2));
      RooChebychev cheb4_bkgPDF_out("cheb4_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb4_x0,bkg_cheb4_x1,bkg_cheb4_x2,bkg_cheb4_x3));
      RooChebychev cheb5_bkgPDF_out("cheb5_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb5_x0,bkg_cheb5_x1,bkg_cheb5_x2,bkg_cheb5_x3,bkg_cheb5_x4));
      RooGamma gamma_bkgPDF_out("gamma_bkgPDF_"+varname,"", dilep_mass_out, bkg_gamma_a,bkg_gamma_b, bkg_gamma_mu);
      RooExponential polexp_bkgPDF1_out("polexp_bkgPDF1_"+varname,"", dilep_mass,bkg_polexp_x0 );
      RooExponential polexp_bkgPDF2_out("polexp_bkgPDF2_"+varname,"", dilep_mass,bkg_polexp_x1 );
     RooAddPdf polexp_bkgPDF_out("polexp_bkgPDF_"+varname,"",RooArgList(polexp_bkgPDF1_out,polexp_bkgPDF2_out),RooArgList(bkg_polexp_c0,bkg_polexp_c1),false);
      RooPolynomial pol4_bkgPDF_out("pol4_bkgPDF_"+varname,"",dilep_mass,RooArgList(bkg_pol4_x0,bkg_pol4_x1, bkg_pol4_x2,bkg_pol4_x3,bkg_pol4_x4,bkg_pol4_x5));
      
      RooCategory cat("pdfindex_"+name, "");
      RooArgList models;
      models.add(cheb3_bkgPDF_out);
      if (altfit_bkg){
        if (altfit_bkg_cheb4) models.add(cheb4_bkgPDF_out);
        if (altfit_bkg_cheb5) models.add(cheb5_bkgPDF_out);
        if (altfit_bkg_bst3) models.add(bst3_bkgPDF_out);
        if (altfit_bkg_bst4) models.add(bst4_bkgPDF_out);
        if (altfit_bkg_gamma) models.add(gamma_bkgPDF_out);
        if (altfit_bkg_exp)  models.add(polexp_bkgPDF_out);
        if (altfit_bkg_pol4)  models.add(pol4_bkgPDF_out);
      }

      RooMultiPdf multipdf("multipdf", "", cat, models);
      RooRealVar norm_out("multipdf_norm","",dataset.numEntries(),0,2*dataset.numEntries());
      
      RooWorkspace *wspace_bkg = new RooWorkspace("workspace_background","workspace_background");
      wspace_bkg->import(cat, RooFit::RecycleConflictNodes());
      wspace_bkg->import(norm_out);
      wspace_bkg->import(multipdf, RooFit::RecycleConflictNodes());
      
      RooDataSet * data_obs = cheb3_bkgPDF_out.generate(RooArgSet(dilep_mass_out),dataset.numEntries());
      cout<<dataset.numEntries()<<endl;
      data_obs->SetName("sim_data_obs_"+varname);
      auto gen_frame = dilep_mass_out.frame();
      data_obs->plotOn(gen_frame,RooFit::Binning(400),RooFit::MarkerColor(1),RooFit::LineColor(1));
      save_plot(gen_frame,"m(e,#mu)","v10_gentest_"+name);
      dilep_mass.setBins(400);
      RooDataHist hdata_sim("sim_binned_obs_"+varname,"sim_binned_obs_"+name,dilep_mass_out,*data_obs);
      hdata_sim.SetName("sim_binned_obs_"+varname);
      RooPlot* bin_frm = dilep_mass_out.frame();
      hdata_sim.plotOn(bin_frm);
      TCanvas * csim = new TCanvas("csin","",700,700);
      bin_frm->Draw();
      csim->SaveAs("v10_datasim_"+name+".png");
      wspace_bkg->import(*data_obs);
      wspace_bkg->import(hdata_sim);
      wspace_bkg->writeToFile("workspace_v10_bkg_"+name+".root");
   }
   
   double nSgn = n_sgn.getVal();
   RooAbsReal* rsgn = sgn_PDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr"));
   double nSgn_inSR = rsgn->getVal()*n_sgn.getVal();

   double nBkg_cheb4 = cheb4_n_bkg.getVal();
   double rBkg_cheb4 = (cheb4_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();

   double nBkg_cheb3 = cheb3_n_bkg.getVal();
   double rBkg_cheb3 = (cheb3_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();

   double nBkg_cheb5 = cheb5_n_bkg.getVal();
   double rBkg_cheb5 = (cheb5_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();

   double nBkg_bst3 = bst3_n_bkg.getVal();
   double rBkg_bst3 = (bst3_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();

   double nBkg_bst4 = bst4_n_bkg.getVal();
   double rBkg_bst4 = (bst4_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();

   double nBkg_gamma = gamma_n_bkg.getVal();
   double rBkg_gamma = (gamma_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();

   double nBkg_polexp = polexp_n_bkg.getVal();
   double rBkg_polexp = (polexp_bkgPDF.createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr")))->getVal();
 
   esgn_PDF.plotOn(xframe,RooFit::LineColor(kBlue),RooFit::Normalization(1, RooAbsReal::RelativeExpected),RooFit::Name("sgn"));
   leg->AddEntry(xframe->findObject("sgn"),"Signal(DCB)");
   save_plot(xframe,"m(e,#mu)","v10_data_fit_and_sgn_"+name,leg);
  
   cout<<"\n\nYieldd for "+name<<endl;


   cout<<"whole range:"<<endl;
   cout<<"Total Evt count: "<<dataset.numEntries()<<endl;
   cout<<" - nSgn "<<nSgn<<endl;
   cout<<" - nBkg(Chebychev3) [main] "<<nBkg_cheb3<<endl;
   if (altfit_bkg_cheb4) cout<<" - nBkg(Chebychev4) "<<nBkg_cheb4<<endl;
   if (altfit_bkg_cheb5) cout<<" - nBkg(Chebychev5) "<<nBkg_cheb5<<endl;
   if (altfit_bkg_bst3) cout<<" - nBkg(Bernstein3) "<<nBkg_bst3<<endl;
   if (altfit_bkg_bst4) cout<<" - nBkg(Bernstein4) "<<nBkg_bst4<<endl;
   if (altfit_bkg_gamma) cout<<" - nBkg(Gamma) "<<nBkg_gamma<<endl;
   if (altfit_bkg_exp) cout<<" - nBkg(Exp Pol) "<<nBkg_polexp<<endl;
   
   
   cout<<"85-95 range:"<<endl;
   cout<<" - nSgn "<<nSgn_inSR<<endl;
   cout<<" - nBkg(Chebychev3) [main] "<<nBkg_cheb3*rBkg_cheb3<<endl;
   if (altfit_bkg_cheb4) cout<<" - nBkg(Chebychev4) "<<nBkg_cheb4*rBkg_cheb4<<endl;
   if (altfit_bkg_cheb5) cout<<" - nBkg(Chebychev5) "<<nBkg_cheb5*rBkg_cheb5<<endl;
   if (altfit_bkg_bst3) cout<<" - nBkg(Bernstein3) "<<nBkg_bst3*rBkg_bst3<<endl;
   if (altfit_bkg_bst4) cout<<" - nBkg(Bernstein4) "<<nBkg_bst4*rBkg_bst4<<endl;
   if (altfit_bkg_gamma) cout<<" - nBkg(Gamma) "<<nBkg_gamma*rBkg_gamma<<endl;
   if (altfit_bkg_exp) cout<<" - nBkg(Exp Pol) "<<nBkg_polexp*rBkg_polexp<<endl;
   
   if (ntoys>0){
     cout<<"pull test"<<endl;
     TH1F * htoy_cheb3_nevt = new TH1F("htoy_cheb3_nevt","",100,dataset.numEntries()-5.*TMath::Sqrt(dataset.numEntries()),dataset.numEntries()+5.*TMath::Sqrt(dataset.numEntries()));
     auto toy_cheb3_datasets = generate_pseudo_data(ntoys, dataset.numEntries(), cheb3_bkgPDF, dilep_mass, true );
     for (int itoy=0; itoy<toy_cheb3_datasets.size(); itoy++){
       htoy_cheb3_nevt->Fill(toy_cheb3_datasets[itoy].numEntries());
     }
     TCanvas * ctoy_nevt = new TCanvas("ctoy_nevt","",700,700);
     htoy_cheb3_nevt->Draw();
     ctoy_nevt->SaveAs("ctoy_nevt.png");
     TH1F * hdiff_bkgonly_exp = new TH1F("hdiff_bkgonly_exp","",100,dataset.numEntries()*0.8,dataset.numEntries()*1.2);
     TH1F * hpull_sgn0_exp = new TH1F("hpull_sgn0_exp","",100,-5,5);
     n_sgn.setConstant(false);
     RooAddPdf total_pdf_exp_bkg("total_pdf_exp_bkg","",RooArgList(polexp_bkgPDF,sgn_PDF),RooArgList(polexp_n_bkg,n_sgn));
//     std::vector<RooFitResult*> toy_ex_results;
     for (auto toy_cheb3: toy_cheb3_datasets ){
        if (altfit_bkg && altfit_bkg_exp){
          bkg_polexp_x0.setVal(-1.8357e-05);
          bkg_polexp_x1.setVal(-9.0861e-02);
          bkg_polexp_c0.setVal(1.e3);
          bkg_polexp_c1.setVal(1.e3);
          polexp_ebkgPDF.fitTo(toy_cheb3,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"));
          total_pdf_exp_bkg.fitTo(toy_cheb3,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"));
          hdiff_bkgonly_exp->Fill( polexp_n_bkg.getVal() - toy_cheb3.numEntries() );
          hpull_sgn0_exp->Fill(n_sgn.getVal()/n_sgn.getError());
//          cout<<n_sgn.getVal()<<"  "<<n_sgn.getError()<<endl;
        }  
     }
     TCanvas * cdiff_bkgonly_exp_cheb3 = new TCanvas("cdiff_bkgonly_exp_cheb3","",700,700);
     hdiff_bkgonly_exp->Draw();
     cdiff_bkgonly_exp_cheb3->SaveAs("cdiff_bkgonly_exp_cheb3.png");
     TCanvas * cpull_sgn0_exp_cheb3 = new TCanvas("cpull_sgn0_exp_cheb3","",700,700);
     hpull_sgn0_exp->Draw();
     cpull_sgn0_exp_cheb3->SaveAs("cpull_sgn0_exp_cheb3.png");

   }
   
 return 0;
} 
