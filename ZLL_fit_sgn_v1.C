#include <fstream>
#include "fit_helper.h"
#include "RooCMSShape.h"


int ZLL_fit_sgn_v1(TString name="test", TString lepton="mu", TString data_file="bdt_tree_v2_data_mumu_part*.root", bool save_shape=true, int nbin=250, bool verbose=false){

   

   gROOT->SetBatch(true);
   TString dilep_var_name="mass_ll";
   TString data_cuts = "mass_ll>70 && mass_ll<110 && pt_ll<45 && mt_l1<70 && mt_l2<70 && ht<50 && met<30 && pt_l2>23";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;

   float sgn_cb_a1_min=0, sgn_cb_a1_max= 10.0; 
   float sgn_cb_n1_min= 0.001, sgn_cb_n1_max=10.0;
   float sgn_cb_n2_min= 0.001, sgn_cb_n2_max=10.0;    
   float cb_width_min=1.0;
   TString leplep="mumu";

   if (lepton=="ele"){
      sgn_cb_a1_min=0.00001;
      sgn_cb_n1_min= 0.1;
      sgn_cb_n1_max=200.0;
      sgn_cb_n2_min=0.1;
      sgn_cb_n2_max=200.0;
      cb_width_min=0.01;
      leplep="ee";
   }

   name = leplep+"_"+name;
   if (!verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);  
   
   TTree * data_tree = get_tree("mytree",data_file,data_cuts);
   RooRealVar dilep_mass("mass_ll","m(l,l)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   
   RooWorkspace ws ("ws");
   ws.import(dilep_mass);
   
   /////////////////////////// Total ////////////////////////////////////////
   RooDataSet dataset("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));
   
   dilep_mass.setBins(nbin);
   RooDataHist datah("data_obs","",RooArgSet(dilep_mass),dataset);
   
   ws.import(dataset); 
   ws.import(datah);
   RooRealVar df("df","",1.0, 0, 1.0);
   df.setConstant(true);  
 
   RooRealVar sgn_cb_mean("sgn_cb_mean","",91.0e+00, 84.0e+00, 96.0e+00);
   RooRealVar sgn_cb_width("sgn_cb_width","",5., cb_width_min, 10.);
   RooRealVar sgn_cb_a1("sgn_cb_a1","",1.0, sgn_cb_a1_min, sgn_cb_a1_max); //for ee 0.0001
   RooRealVar sgn_cb_n1("sgn_cb_n1","",1.0, sgn_cb_n1_min, sgn_cb_n1_max); //for 0.1, 200.0);
   RooRealVar sgn_cb_a2("sgn_cb_a2","",1.0, 0.0001, 10.0); //for ee 
   RooRealVar sgn_cb_n2("sgn_cb_n2","",1.0, sgn_cb_n2_min, sgn_cb_n2_max); // for ee 1.0, 0.1, 200.0);
   RooCrystalBall sgn_cb1("sgn_cb1","cb",dilep_mass,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);

   RooRealVar sgn_g_mean("sgn_g_mean","",91.0e+00, 84.0e+00, 96.0e+00);
   RooRealVar sgn_g_width("sgn_g_width","",5., 0.1, 10.);
   RooGaussian sgn_g("sgn_g","sgn_g",dilep_mass,sgn_g_mean,sgn_g_width);

   RooRealVar sgn_frac("sgn_frac","",0.5, 0, 10.0);
   RooAddPdf sgnPDF("sgnPDF","",RooArgList(sgn_cb1, sgn_g),RooArgList(df,sgn_frac));
   RooRealVar nsgn("nsgn","",10000,0,1000000000000);
   RooExtendPdf total_PDF("total_PDF","",sgnPDF,nsgn);

   RooFitResult * total_result = total_PDF.fitTo(datah,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1));
   
   auto total_frame = dilep_mass.frame();
   datah.plotOn(total_frame);
   total_PDF.plotOn(total_frame,RooFit::LineColor(kRed),RooFit::Normalization(1, RooAbsReal::RelativeExpected));
   print_details (total_result, total_frame);
   save_plot(total_frame,"m(l,l)","v1sgn_cr_"+name,true);
   save_pull(total_frame, dilep_mass, "m(l,l)","v1sgn_cr_"+name);
   cout<<"Sgn "<<nsgn.getVal()<<endl;

   if (save_shape){
      TFile * fHist = new TFile("cr_v1fitsgn_shape_"+name+"-TH1.root", "RECREATE");
      RooDataSet * signal_sim = sgnPDF.generate(RooArgSet(dilep_mass),nsgn.getVal());
      TH1 * histo =  sgnPDF.createHistogram(leplep,dilep_mass);
      histo->SetName(leplep);
//      histo->Scale(nsgn.getVal()/histo->Integral());
      TH1 * hdata_out = datah.createHistogram("data_obs",dilep_mass);
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
