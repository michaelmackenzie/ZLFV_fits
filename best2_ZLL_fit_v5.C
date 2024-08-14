#include <fstream>
#include "fit_helper.h"
#include "RooCMSShape.h"
R__LOAD_LIBRARY(RooCMSShape_cc.so)




int ZLL_fit_v5(TString name="test", TString lepton="mu", TString data_file="bdt_tree_v2_data_mumu_part*.root", bool save_shape=true, int nbin=500, bool save_data_histo=true, TString directly_input_histo_file="None",TString input_histo="None"){
   
   TString dilep_var_name="mass_ll";
   TString cuts = "mass_ll>70 && mass_ll<110 && xgb>0.65";
   TString tree_name="mytreefit";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;
   bool verbose=false;   




   TString leplep="mumu";
   if (lepton=="ele"){
      leplep="ee";
   }
   
   if (!verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);  

   TH1F* hdata = new TH1F("hdata","",nbin,min_fit_range,max_fit_range);
   if (directly_input_histo_file=="None"){
      cout<<"Histo is created from "+data_file+" with defined cuts/bins"<<endl;
      TTree * data_tree = get_tree(tree_name,data_file,cuts);
      data_tree->Draw(dilep_var_name+">>hdata",cuts);
      if (save_data_histo){
        TFile * fdata_histo = new TFile("fdata_from_fit.root","RECREATE");
        hdata->Write();
        cout<<"====> saved in fdata_from_fit.root"<<endl;
      }
   } else{
     cout<<" Reading DIRECTLY stored histo ("+input_histo+") from file "+directly_input_histo_file+"; cuts/data-file/bins NOT used"<<endl;  
     TFile * fin = new TFile(directly_input_histo_file,"READ");
     hdata = (TH1F*) fin->Get(input_histo);
//     hdata->Rebin(4);
   }

   RooRealVar dilep_mass("mass_ll","m(l,l)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   
   RooWorkspace ws ("ws");
   ws.import(dilep_mass);
   RooRealVar df("df","",1.0, 0, 1.0);
   df.setConstant(true);

//   RooDataSet dataset("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));

   
   dilep_mass.setBins(nbin);
//   RooDataHist datah("datah","",RooArgSet(dilep_mass),dataset);
   RooDataHist datah("datah","",RooArgSet(dilep_mass),hdata);

   RooRealVar bkg_cmsshape_a("bkg_cmsshape_a", "", 60., 50., 190.); // 60., 50., 190.
   RooRealVar bkg_cmsshape_b("bkg_cmsshape_b", "", 0.05, 0.01, 0.24); // 0.05, 0.01, 0.12
   RooRealVar bkg_cmsshape_c("bkg_cmsshape_c", "", 0.1, -2., 2.);
   RooRealVar bkg_cmsshape_peak("bkg_cmsshape_peak", "", 90); 

   RooCMSShape bkg_cmsshape("bkg_cmsshape", "", dilep_mass,bkg_cmsshape_a,bkg_cmsshape_b,bkg_cmsshape_c,bkg_cmsshape_peak);

   RooRealVar sgn_voi1_mean("sgn_voi1_mean", "", 91.0,90.,92.); //90
   RooRealVar sgn_voi1_width("sgn_voi1_width", "", 2.495);
   RooRealVar sgn_voi1_sigma("sgn_voi1_sigma", "", 0.9, 0.5, 5.); //0.8 0.05
   //RooVoigtian sgn_voi1("sgn_voi1","",dilep_mass,sgn_voi1_mean,sgn_voi1_width,sgn_voi1_sigma);
   RooVoigtian sgn_voi1("sgn_voi1","",dilep_mass,sgn_voi1_mean,sgn_voi1_width,sgn_voi1_sigma);

   RooRealVar sgn_cb_mean("sgn_cb_mean", "", 0., -5., 5.);
   RooRealVar sgn_cb_sigma("sgn_cb_sigma", "", 0.9,0.005,5.0);
   RooRealVar sgn_cb_a("sgn_cb_a", "", 4.0, 0.0, 10.);
   RooRealVar sgn_cb_n("sgn_cb_n", "", 4.0, 0.0, 10.);
   RooCBShape sgn_cb("sgn_cb","",dilep_mass,sgn_cb_mean,sgn_cb_sigma,sgn_cb_a,sgn_cb_n);

//   RooRealVar sgn_frac("sgn_frac","",0.7,0.01,10.0);
//   RooAddPdf sgn_voi1("sgn_voi1","",RooArgList(sgn_voi1tmp, sgn_erf),RooArgList(df,sgn_frac));

   dilep_mass.setBins(10000,"cache") ;
   RooFFTConvPdf sgn_PDF("sgn_PDF","",dilep_mass,sgn_voi1,sgn_cb) ;

   RooRealVar nbkg("nbkg","",10000,0,1000000000000);
   RooRealVar nsgn("nsgn","",10000,0,1000000000000);


   RooAddPdf total_PDF("total_PDF","",RooArgList(sgn_PDF, bkg_cmsshape),RooArgList(nsgn,nbkg));
   
   RooFitResult * result = total_PDF.fitTo(datah,RooFit::Extended(1),RooFit::Save(),RooFit::Range("fitRange"),RooFit::Minos(false)); //RooFit::PrintLevel(-1)

   auto frame = dilep_mass.frame();
   datah.plotOn(frame);
   sgn_PDF.plotOn(frame,RooFit::LineColor(kBlue),RooFit::LineStyle(1),RooFit::Normalization(nsgn.getVal(), RooAbsReal::NumEvent));
   bkg_cmsshape.plotOn(frame,RooFit::LineColor(kBlack),RooFit::LineStyle(1),RooFit::Normalization(nbkg.getVal(), RooAbsReal::NumEvent));
   total_PDF.plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(1, RooAbsReal::RelativeExpected));
   print_details (result, frame);
   save_plot(frame,"m(l,l)","v5_cr_"+leplep+"_"+name,true,false);
   save_pull(frame, dilep_mass, "m(l,l)","v4_cr_"+leplep+"_"+name);
   cout<<endl<<"RESULTS:"<<endl;
//   cout<<"  -- count evts "<<datah.GetEntries()<<endl;
   cout<<"  -- Bkg "<<nbkg.getVal()<<endl;
   cout<<"  -- Sgn "<<nsgn.getVal()<<endl<<endl;

   if (save_shape){
      TFile * fHist = new TFile("cr_v5fit_shape_"+leplep+"_"+name+"-TH1.root", "RECREATE");
      RooDataSet * signal_sim = sgn_PDF.generate(RooArgSet(dilep_mass),nsgn.getVal());
      RooDataSet * data_sim = sgn_PDF.generate(RooArgSet(dilep_mass),nsgn.getVal());
      TH1 * histo =  signal_sim->createHistogram(leplep,dilep_mass,nbin);
      histo->SetName(leplep);
      TH1 * hdata_out = data_sim->createHistogram("data_obs",dilep_mass,nbin);
      TCanvas * csgn_check = new TCanvas("csgn_check","",700,700);
      histo->Scale(hdata_out->Integral()/histo->Integral());
      histo->Draw("HIST");
      histo->Write();
      cout<<"Sgn2 "<<hdata_out->Integral()<<endl;
      hdata_out->SetName("data_obs");
      hdata_out->Write();
   }

 return 0;
} 
