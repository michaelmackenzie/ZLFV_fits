#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// mk2: seperating fit in components
///////////// Continuous bkg component v1 
//////////////////// + Same thing as v401 etc
//////////////////// + Added Gauss&Expo
//////////////////// + Added Gauss&Power
//////////////////// + Ability to write families in ws



int ZMuE_fit_mk2_bkg_v1(TString name="bin1_r2", 
    TString bkg_file="pseudo_data_from_MC_v2_r0.root",
    TString xgbmin="0.3",TString xgbmax="0.7", bool create_dc_input=false, 
    TString outvar="mass_ll", TString varname="bin", 
    bool pseudodata_input=false, float pseudodata_norm=-1.0, 
    bool histo_input=false){

   //////////////////////////////////// configuration /////////////////////////
   gROOT->SetBatch(true);
   TString ztt_file = "/afs/cern.ch/work/g/gkaratha/private/Analysis/DispJets/Analyzer/CMSSW_10_2_16_UL/src/PhysicsTools/NanoAODTools/plotter/ZMuE_plotting_and_cfg/BDT/Meas_fullAndSFAndGenDecay_bdt_v7_bkg_dy_mcRun2_extend.root";


   int min_gspol_order=1,max_gspol_order=3;
   bool Fit_gspol=true;
   int add_orders_gspol = 0;
   int min_gsexp_order=1,max_gsexp_order=3;
   bool Fit_gsexp=true;
   int add_orders_gsexp = 0;
   int min_gsplaw_order=1,max_gsplaw_order=3;
   bool Fit_gsplaw=true;
   int add_orders_gsplaw = 0;
 
   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron"; // not add mass_ll here when run systematics
   TString tree_name="mytreefit";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;
   double blind_min = 86.;
   double blind_max = 96.;
   int nbin_data = 80;
   bool Verbose=false; // RooFit verbosity
   bool Bkg_only_fit_whole_region=false;
   int printout_levels=1; // 0: Print only final fits parameters, 1: Print all fits from F test tests

   

   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ZMuE mk2 fit v1: bkg part ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;
   
   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   

   // read trees
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   dilep_mass.setBins(nbin_data);
   dilep_mass_out.setBins(nbin_data);

   dilep_mass.setRange("left",min_fit_range, blind_min); //70,88
   dilep_mass.setRange("right",blind_max, max_fit_range); //94,110
   dilep_mass.setRange("sr",blind_min,blind_max); //88,94
   dilep_mass.setRange("full",min_fit_range, max_fit_range);


   ///////////////////////////////// Data/bkg fit ////////////////////////////
   cout<<"\n *********************** Z->tt bump fit ********************** "<<endl;

   //////// First fit on the Z->tt bump bellow 70 GeV -> gaussian part
   TChain * ccdy = new TChain("mytreefit");
   ccdy->Add(ztt_file);

   TH1F * hmcext_ztt = new TH1F("hmcext_ztt","",100,50,150);
   ccdy->Draw("mass_ll>>hmcext_ztt",xgbmin+"<xgb && xgb<"+xgbmax+" && IsGen_Ztt==1");
 
   RooRealVar mll_mc_extd("mass_ll","m(e,#mu)", 70., 50, 80, "GeV/c^{2}");
   mll_mc_extd.setBins(30);
   RooDataHist dhistext_ztt("dhistext_ztt","",RooArgSet(mll_mc_extd),hmcext_ztt);
   
   RooRealVar * ztt_gs_mu= new RooRealVar("ztt_gs_mu_"+varname,"ztt_gs_mu_"+varname,60,50,69); //69
   RooRealVar * ztt_gs_wd = new RooRealVar("ztt_gs_wd_"+varname,"ztt_gs_wd_"+varname,11,5,20);

   RooGaussian *ztt_gs_pdf =new RooGaussian("ztt_gs_pdf","ztt_gs_pdf",mll_mc_extd,*ztt_gs_mu,*ztt_gs_wd);

   RooRealVar ztt_gs_pdf_norm("ztt_gs_pdf_norm", "ztt_gs_pdf_norm",0,0,2*dhistext_ztt.sumEntries());
   RooAddPdf ztt_gs_epdf("ztt_gs_epdf","",RooArgList(*ztt_gs_pdf),RooArgList(ztt_gs_pdf_norm));

   ///// fit
   RooFitResult * ztt_gs_result = ztt_gs_epdf.fitTo(dhistext_ztt,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1));
   
   print_details (ztt_gs_result);

   auto ext_ztt_frame = mll_mc_extd.frame();
   dhistext_ztt.plotOn(ext_ztt_frame,RooFit::Binning(30),RooFit::Name("data"),RooFit::LineColor(1),RooFit::MarkerColor(1));
   ztt_gs_pdf->plotOn(ext_ztt_frame,RooFit::Binning(30),RooFit::Name("ztt_gs_pdf"),RooFit::LineColor(1),RooFit::MarkerColor(1));
   save_plot(ext_ztt_frame,"m(e,#mu)","v401_ztt_gauss_"+name,new TLegend(),NULL, true, true);
   save_plot_and_band(ext_ztt_frame,mll_mc_extd,{"ztt_gs_pdf"},"m(e,#mu)","v401_ztt_gauss_band_"+name);


    
   ////// get bkg dataset
   cout<<"\n *********************** Background fit ********************** "<<endl;

   TH1F * hbkg = new TH1F("hbkg","",nbin_data,min_fit_range,max_fit_range);

   if (histo_input){
     TFile * fzmm = new TFile(bkg_file,"READ");
     THStack * stzmm = (THStack *) fzmm->Get("bkg_stack");
     auto * hists = stzmm -> GetHists();
     for (auto * hist: *hists){
       TString hist_name =hist->GetName();
       TH1F* th1 = (TH1F*) hist;
       if (!hist_name.BeginsWith("Z"))
          hbkg->Add(th1);
     }
   } else {
     TChain * cc = new TChain(tree_name);
     cc->Add(bkg_file);
     if (pseudodata_norm>0) 
       cc->Draw("mass_ll>>hbkg",TString(std::to_string(pseudodata_norm))+"*NormGen_wt*( !IsZmm  && !IsZee && "+cuts+" )");
     else
       cc->Draw("mass_ll>>hbkg","NormGen_wt*( !IsZmm  && !IsZee && "+cuts+" )");
   }

   RooDataHist * dhist_bkg = new RooDataHist("dhist_bkg","dhist_bkg",RooArgSet(dilep_mass),hbkg);


   /////// gaus + polynomial 
   std::vector<RooAbsPdf*> bkg_gs_pol_pdfs;
   std::vector<RooRealVar*> bkg_gs_pol_ampl;
   std::vector<TString> bkg_gs_pol_names;
   std::vector<TString> bkg_gs_pol_legs;
   std::vector<int> bkg_gs_pol_orders;

   for (int iorder=min_gspol_order; iorder<max_gspol_order+1; iorder++){
      TString sorder(std::to_string(iorder));
      bkg_gs_pol_pdfs.push_back(CreateGaussPolynomial( "gauss_pol"+sorder+"_pdf", iorder, dilep_mass,ztt_gs_mu,ztt_gs_wd));
      bkg_gs_pol_ampl.push_back(new RooRealVar("gauss_pol"+sorder+"_pdf_norm", "gauss_pol"+sorder+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
      bkg_gs_pol_names.push_back("gauss_pol"+sorder+"_"+name);
      bkg_gs_pol_legs.push_back("Gauss+Polynomial "+sorder);
      bkg_gs_pol_orders.push_back(iorder);
   }

   ///// ftest
   FtestStruct gs_pol_Ftest;
   gs_pol_Ftest.success=false;
   if (Fit_gspol) {
      cout<<" ************************ No Z->mumu Ftest Gauss + Poly begin ************************ "<<endl;
      gs_pol_Ftest =  HistFtest(bkg_gs_pol_pdfs, bkg_gs_pol_ampl,  dhist_bkg, dilep_mass, bkg_gs_pol_orders, bkg_gs_pol_names, bkg_gs_pol_legs, nbin_data,"mk2bkg_v1_gauspol_"+name, 0.05,0.01,printout_levels);
      cout<<" ************************ No Z->mumu Ftest Gauss + Poly end ************************ "<<endl;
   }

   /////// gaus + expo
   std::vector<RooAbsPdf*> bkg_gs_exp_pdfs;
   std::vector<RooRealVar*> bkg_gs_exp_ampl;
   std::vector<TString> bkg_gs_exp_names;
   std::vector<TString> bkg_gs_exp_legs;
   std::vector<int> bkg_gs_exp_orders;

   for (int iorder=min_gsexp_order; iorder<max_gsexp_order+1; iorder++){
      TString sorder(std::to_string(iorder));
      bkg_gs_exp_pdfs.push_back(CreateGaussExpo( "gauss_exp"+sorder+"_pdf", iorder, dilep_mass,ztt_gs_mu,ztt_gs_wd));
      bkg_gs_exp_ampl.push_back(new RooRealVar("gauss_exp"+sorder+"_pdf_norm", "gauss_exp"+sorder+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
      bkg_gs_exp_names.push_back("gauss_exp"+sorder+"_"+name);
      bkg_gs_exp_legs.push_back("Gauss+Expo "+sorder);
      bkg_gs_exp_orders.push_back(iorder);
   }

   ///// ftest
   FtestStruct gs_exp_Ftest;
   if (Fit_gsexp) {
      cout<<" ************************ No Z->mumu Ftest Gauss + Expo begin ************************ "<<endl;
      gs_exp_Ftest =  HistFtest(bkg_gs_exp_pdfs, bkg_gs_exp_ampl,  dhist_bkg, dilep_mass, bkg_gs_exp_orders, bkg_gs_exp_names, bkg_gs_exp_legs, nbin_data, "mk2bkg_v1_gausexp_"+name, 0.05,0.01,printout_levels);
      cout<<" ************************ No Z->mumu Ftest Gauss + Expo end ************************ "<<endl;
   }

   /////// gaus + power law
   std::vector<RooAbsPdf*> bkg_gs_plaw_pdfs;
   std::vector<RooRealVar*> bkg_gs_plaw_ampl;
   std::vector<TString> bkg_gs_plaw_names;
   std::vector<TString> bkg_gs_plaw_legs;
   std::vector<int> bkg_gs_plaw_orders;

   for (int iorder=min_gsplaw_order; iorder<max_gsplaw_order+1; iorder++){
      TString sorder(std::to_string(iorder));
      bkg_gs_plaw_pdfs.push_back(CreateGaussPower( "gauss_plaw"+sorder+"_pdf", iorder, dilep_mass,ztt_gs_mu,ztt_gs_wd));
      bkg_gs_plaw_ampl.push_back(new RooRealVar("gauss_plaw"+sorder+"_pdf_norm", "gauss_plaw"+sorder+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
      bkg_gs_plaw_names.push_back("gauss_plaw"+sorder+"_"+name);
      bkg_gs_plaw_legs.push_back("Gauss+Power law "+sorder);
      bkg_gs_plaw_orders.push_back(iorder);
   }
   ///// ftest
   FtestStruct gs_plaw_Ftest;
   if (Fit_gsplaw){
     cout<<" ************************ No Z->mumu Ftest Gauss + Power law begin ************************ "<<endl;
     gs_plaw_Ftest =  HistFtest(bkg_gs_plaw_pdfs, bkg_gs_plaw_ampl,  dhist_bkg, dilep_mass, bkg_gs_plaw_orders, bkg_gs_plaw_names, bkg_gs_plaw_legs, nbin_data, "mk2bkg_v1_gausplaw_"+name, 0.05,0.01,printout_levels);
     cout<<" ************************ No Z->mumu Ftest Gauss + Power law end ************************ "<<endl;
   }
   
   //////////////////// fit all candidate functions ////////////////
   std::vector<RooAbsPdf*> bkg_pdfs;
   std::vector<RooRealVar*> bkg_ampl;
   std::vector<TString> bkg_names;
   std::vector<TString> bkg_legs;

   if (gs_pol_Ftest.success) {
     for (int i=0; i<gs_pol_Ftest.getAllOrder.size(); i++){
        int iord =gs_pol_Ftest.getAllOrder[i];
        TString sord(std::to_string(iord));
        bkg_pdfs.push_back(CreateGaussPolynomial( varname+"_gspol"+sord+"_pdf", iord, dilep_mass,ztt_gs_mu,ztt_gs_wd));
        bkg_ampl.push_back(new RooRealVar(varname+"_gspol"+sord+"_pdf_norm", varname+"_gspol"+sord+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
        bkg_names.push_back(varname+"_gspol"+sord+"_pdf");
        bkg_legs.push_back("Gauss+Polynomial "+sord);
     }
   }
   
   if (gs_exp_Ftest.success) {
     for (int i=0; i<gs_exp_Ftest.getAllOrder.size(); i++){
        int iord =gs_exp_Ftest.getAllOrder[i];
        TString sord(std::to_string(iord));
        bkg_pdfs.push_back(CreateGaussExpo( varname+"_gsexp"+sord+"_pdf", iord, dilep_mass,ztt_gs_mu,ztt_gs_wd));
        bkg_ampl.push_back(new RooRealVar(varname+"_gsexp"+sord+"_pdf_norm", varname+"_gsexp"+sord+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
        bkg_names.push_back(varname+"_gsexp"+sord+"_pdf");
        bkg_legs.push_back("Gauss+Expo "+sord);
     }
   }
   
   if (gs_plaw_Ftest.success) {
     for (int i=0; i<gs_plaw_Ftest.getAllOrder.size(); i++){
        int iord =gs_plaw_Ftest.getAllOrder[i];
        TString sord(std::to_string(iord));
        bkg_pdfs.push_back(CreateGaussPower( varname+"_gsplaw"+sord+"_pdf", iord, dilep_mass,ztt_gs_mu,ztt_gs_wd));
        bkg_ampl.push_back(new RooRealVar(varname+"_gsplaw"+sord+"_pdf_norm", varname+"_gsplaw"+sord+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
        bkg_names.push_back(varname+"_gsplaw"+sord+"_pdf");
        bkg_legs.push_back("Gauss+Power "+sord);
     }
   }



   //////////// fit here
   cout<<" Fit of best variables "<<endl;
   std::vector<std::vector<float>> final_results =  FitHistBkgFunctions(bkg_pdfs, bkg_ampl, dhist_bkg, dilep_mass, bkg_names, bkg_legs, nbin_data, true, "mk2bkg_v1_best_"+name, true);

   //////////////////// create bkg workspace ////////////////////////
   if (!create_dc_input)
   return 0;

   RooCategory cat("pdfindex_"+varname, "");
   RooArgList models_out; //container for multpdf
   RooWorkspace *wspace = new RooWorkspace("workspace_background","workspace_background"); //background & data workspace

   int nstartFNC1=0;
   if (gs_pol_Ftest.success) {
     for (int i=0; i<gs_pol_Ftest.getAllOrder.size(); i++){
       std::vector<float> param;
       TString sord (std::to_string(gs_pol_Ftest.getAllOrder[i]));
       for(int j=2; j<final_results[i].size(); j++)
          param.push_back(final_results[i][j]);
       models_out.add( *(CreateGaussPolynomial( "bkg_gspol"+sord+"_pdf_"+varname, gs_pol_Ftest.getAllOrder[i], dilep_mass_out,ztt_gs_mu,ztt_gs_wd, param)) );
       nstartFNC1+=1;
     }
     for( int iord =gs_pol_Ftest.getBestOrder+1; iord<gs_pol_Ftest.getBestOrder+add_orders_gspol+1; iord++){
       TString sord (std::to_string(iord));
       models_out.add( *(CreateGaussPolynomial( "bkg_gspol"+sord+"_pdf_"+varname, iord, dilep_mass_out,ztt_gs_mu,ztt_gs_wd)) );
     }
   }

   int nstartFNC2=nstartFNC1;
   if (gs_exp_Ftest.success) {
     for (int i=0; i<gs_exp_Ftest.getAllOrder.size(); i++){   
       std::vector<float> param;
       for(int j=2; j<final_results[i+nstartFNC1].size(); j++)
          param.push_back(final_results[i+nstartFNC1][j]);
       TString sord (std::to_string(gs_exp_Ftest.getAllOrder[i]));
       models_out.add( *(CreateGaussExpo( "bkg_gsexp"+sord+"_pdf_"+varname, gs_exp_Ftest.getAllOrder[i], dilep_mass_out,ztt_gs_mu,ztt_gs_wd,param)) );
       nstartFNC2+=1;   
     }  
     for( int iord =gs_exp_Ftest.getBestOrder+1; iord<gs_exp_Ftest.getBestOrder+add_orders_gsexp+1; iord++){
       TString sord (std::to_string(iord));
       models_out.add( *(CreateGaussExpo( "bkg_gsexp"+sord+"_pdf_"+varname, iord, dilep_mass_out,ztt_gs_mu,ztt_gs_wd)) );
     }
   }

   if (gs_plaw_Ftest.success) {
     for (int i=0; i<gs_plaw_Ftest.getAllOrder.size(); i++){   
       std::vector<float> param;
       for(int j=2; j<final_results[i+nstartFNC2].size(); j++)
          param.push_back(final_results[i+nstartFNC2][j]);
       TString sord (std::to_string(gs_plaw_Ftest.getAllOrder[i]));
       models_out.add( *(CreateGaussPower( "bkg_gsplaw"+sord+"_pdf_"+varname, gs_plaw_Ftest.getAllOrder[i], dilep_mass_out,ztt_gs_mu,ztt_gs_wd,param)) );
     }  
    for( int iord =gs_plaw_Ftest.getBestOrder+1; iord<gs_plaw_Ftest.getBestOrder+add_orders_gsplaw+1; iord++){
       TString sord (std::to_string(iord));
       models_out.add( *(CreateGaussPower( "bkg_gsplaw"+sord+"_pdf_"+varname, iord, dilep_mass_out,ztt_gs_mu,ztt_gs_wd)) );
     }
   }
   RooMultiPdf multipdf("multipdf_"+varname, "", cat, models_out);
   RooRealVar norm_out("multipdf_"+varname+"_norm","",dhist_bkg->sumEntries(),0.5*dhist_bkg->sumEntries(),10*dhist_bkg->sumEntries());  
    
   wspace->import(cat);
   wspace->import(multipdf);
   wspace->import(norm_out);
   wspace->writeToFile("workspace_mk2bkg_v1_"+name+".root"); // write outputt   

 return 0;
} 
