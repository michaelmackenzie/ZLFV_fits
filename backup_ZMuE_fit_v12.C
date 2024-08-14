#include <fstream>
#include "fit_helper.h"

//////////  set multiple values in parameter vector
void DefaultVals(std::vector<RooRealVar> vars, std::vector<double> vals){
  for (int i =0; i< vars.size(); i++)
    vars[i].setVal(vals[i]);
}

/////// calculates yield in a subregion and inclusive
std::pair<double,double> yield_calc( float nYld_total, RooRealVar dilep_mass, RooAbsPdf *pdf){
   RooAbsReal* rsgn = pdf->createIntegral(dilep_mass,RooFit::NormSet(dilep_mass),RooFit::Range("sr"));
   double nYld_inSR = rsgn->getVal()*nYld_total;
   return make_pair(nYld_total,nYld_inSR);
}


///////// systematic functions parameters : returns maximum deviation of mean/width
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
    save_plot(syst_frame,"m(#mu,e)","v12_prmtv_syst_"+syst_name);
    delete syst_result;
    return std::make_pair(mean,width);
}


/////////// create cheby parameters
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

////////// create bernstein function parameters
std::vector<RooRealVar> BstParams(int order, TString varname){
   std::vector<RooRealVar> bkg_bsts;
   for (int i =0; i<order; i++){
     RooRealVar bkg_bst_x("bkg_bst"+TString(std::to_string(order))+"_x"+TString(std::to_string(i))+"_"+varname, "", 1./TMath::Power(10,order), -25, 25);
     bkg_bsts.push_back(bkg_bst_x);
   }
   return  bkg_bsts;
}

void PlotFunctions(std::vector<RooAbsPdf*> pdfs, RooPlot * xframe, RooRealVar dilep_mass, std::vector<TString> names, std::vector<TString> legs, TLegend *leg, RooDataSet dataset, TString name, int nbin_data, bool unblind){
   TString norm_range="left,right";
   if (unblind)
      norm_range="full";
   for (int i=0; i<pdfs.size(); i++){
     pdfs[i]->plotOn(xframe,RooFit::LineColor(i+1),RooFit::Range("full"),RooFit::NormRange(norm_range),RooFit::Name(names[i]));
     save_pull(xframe, dilep_mass, "m(e,#mu)","v12_"+names[i]);
     leg->AddEntry(xframe->findObject(names[i]),legs[i]);
   }
   if (unblind)
     dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::Name("data"));
   else
     dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::CutRange("left,right"),RooFit::Name("data"));
   leg->AddEntry(xframe->findObject("data"),"Data");
   save_plot(xframe,"m(e,#mu)","v12_sr_total_"+name,leg);
}

////////// Fit bokg only functions sideband or total fit and blind/unblind plot
void FitBkgFunctions(std::vector<RooAbsPdf*> pdfs, std::vector<RooRealVar*> ampls, RooDataSet &dataset, RooRealVar dilep_mass, vector<TString> names, std::vector<TString> legs, bool allrange, bool unblind, int nbin_data, bool is_pseudodata , TString extra_name){
  TString fit_range="left,right";
  if (allrange) 
     fit_range="full";
  TLegend * leg = new TLegend();
  auto plot_frame =  dilep_mass.frame();
  if (unblind)
    dataset.plotOn(plot_frame,RooFit::Binning(nbin_data),RooFit::Name("data"),RooFit::AsymptoticError(is_pseudodata));
  else  
    dataset.plotOn(plot_frame,RooFit::Binning(nbin_data),RooFit::CutRange("left,right"),RooFit::Name("data"),RooFit::AsymptoticError(is_pseudodata) );

  TPaveText * pt = new TPaveText(0.4,0.8,0.9,0.9,"tlNDC");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->SetTextAlign(12);
  pt->AddText("#chi^{2}:\n");
  cout<<"\n *********** Bkg-only fit in "+fit_range<<" *********"<<endl;
  
  for (int i=0; i<pdfs.size(); i++){ 
    cout<<" Fit: "+names[i]<<endl;
    RooAddPdf epdf("epdf_"+names[i],"", RooArgList(*pdfs[i]),  RooArgList(*ampls[i]));
    RooFitResult * fit_result = epdf.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range(fit_range),RooFit::AsymptoticError(is_pseudodata));
    auto tmp_frame = dilep_mass.frame();
    dataset.plotOn(tmp_frame,RooFit::Binning(nbin_data));
    pdfs[i]->plotOn(tmp_frame,RooFit::LineColor(i+1),RooFit::Range("full"),RooFit::NormRange("full"),RooFit::Name(names[i]));
    print_details (fit_result, tmp_frame);
    int n_param = fit_result->floatParsFinal().getSize();
    float chi2=int(tmp_frame->chiSquare(n_param-2)*100)/100.;
    pdfs[i]->plotOn(plot_frame,RooFit::LineColor(i+1),RooFit::Range("full"),RooFit::NormRange(fit_range),RooFit::Name(names[i]));
    leg->AddEntry(tmp_frame->findObject(names[i]),legs[i]);
    pt->AddText(legs[i]+" "+TString(to_string(chi2)));
    cout<<" ******** Bkg-only result ********* "<<endl;
    cout<<"whole range    |    85-95 only"<<endl;
    std::pair<double,double> nBkg = yield_calc( ampls[i]->getVal(), dilep_mass, pdfs[i]);
    cout<<" - nBkg("+names[i]+") "<<nBkg.first<<"  |  "<<nBkg.second<<endl;
     
  }
  save_plot(plot_frame,"m(e,#mu)","v12_bkgfit_"+extra_name,leg,pt);
  cout<<" ***************************************"<<endl;
}

/////////// adds bkg and signal functions and fit them unblind
void FitTotalFunctions(std::vector<RooAbsPdf*> bkg_pdfs, RooAbsPdf* signal_pdf, RooRealVar dilep_mass, RooDataSet &dataset, vector<TString> names, int nbin_data, bool is_pseudodata, TString extra_name){
  cout<<"\n ************ Unblind fit ************** "<<endl;
  for (int i=0; i<bkg_pdfs.size(); i++){
    cout<<" PDF: "<<names[i]<<endl;
    RooRealVar nBkg_tmp("nBkg_"+names[i],"", dataset.numEntries(),-2*dataset.numEntries() , 2*dataset.numEntries());
    RooRealVar nSgn_tmp("nSgn_"+names[i],"", 0, -dataset.numEntries()/10 , dataset.numEntries()/10);
    RooAddPdf total_pdf("total_pdf_"+names[i],"", RooArgList(*bkg_pdfs[i], *signal_pdf),  RooArgList(nBkg_tmp,nSgn_tmp));
    RooFitResult * fit_result = total_pdf.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"),RooFit::AsymptoticError(is_pseudodata));
    //,RooFit::SumW2Error(true));
    //,RooFit::AsymptoticError(is_pseudodata));
    auto tmpframe = dilep_mass.frame();
    dataset.plotOn(tmpframe,RooFit::Binning(nbin_data),RooFit::AsymptoticError(is_pseudodata));   
    total_pdf.plotOn(tmpframe,RooFit::LineColor(kBlack),RooFit::Range("full"),RooFit::NormRange("full"),RooFit::Name("total_pdf_"+names[i]));
    print_details (fit_result, tmpframe);
    int n_param = fit_result->floatParsFinal().getSize();
    float chi2=int(tmpframe->chiSquare(n_param-2)*100)/100.;
    bkg_pdfs[i]->plotOn(tmpframe,RooFit::LineColor(kRed),RooFit::Normalization(nBkg_tmp.getVal(), RooAbsReal::NumEvent),RooFit::Name("background_pdf_"+names[i]),RooFit::LineStyle(7));
    signal_pdf->plotOn(tmpframe,RooFit::LineColor(kBlue),RooFit::Normalization(nSgn_tmp.getVal(), RooAbsReal::NumEvent),RooFit::Name("signal_pdf_"+names[i]),RooFit::LineStyle(7));
    cout<<" ******** UNBLIND result ********* "<<endl;
    cout<<" - nSgn "<<nSgn_tmp.getVal()<<" +/- "<<nSgn_tmp.getError()<<endl;
    cout<<" - nBkg "<<nBkg_tmp.getVal()<<" +/- "<<nBkg_tmp.getError()<<endl;
    cout<<" ****************************** "<<endl;
    TLegend * leg = new TLegend(0.6,0.6,0.9,0.85);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry(tmpframe->findObject("total_pdf_"+names[i]),"Total");
    leg->AddEntry(tmpframe->findObject("signal_pdf_"+names[i]),"Signal");
    leg->AddEntry(tmpframe->findObject("background_pdf_"+names[i]),"Background");
    TPaveText * pt = new TPaveText(0.6,0.4,0.9,0.6,"brNDC");
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    pt->SetTextAlign(12);
    pt->AddText("SGN yield: "+TString(std::to_string( int(nSgn_tmp.getVal()) ))+" #pm"+TString(std::to_string( int(nSgn_tmp.getError()) )));
    pt->AddText("BKG yield: "+TString(std::to_string( int(nBkg_tmp.getVal()) ))+" #pm"+TString(std::to_string( int(nBkg_tmp.getError()) )));
    pt->AddText("#chi^{2}: "+TString(std::to_string( chi2 )));
    save_plot(tmpframe,"m(e,#mu)","v12_totalfit_"+names[i],leg,pt);
    delete leg;
  }
  cout<<" ****************************** "<<endl;
}


///////// get dataset from data/MC and rescales the latter with weight
RooDataSet GetDataSet(TTree *data_tree, RooRealVar &dilep_mass, bool Is_pseudodata, RooRealVar &wgt, RooFormulaVar &wgtFunc, float norm_factor){
  if (Is_pseudodata){
    TH1F* nbkg = new TH1F("nbkg","",1,0,2); 
    TH1F* nsgn = new TH1F("nsgn","",1,0,2);
    TH1F* ntot = new TH1F("ntot","",1,0,2);
    data_tree->Draw("1>>ntot","NormGen_wt");
    data_tree->Draw("1>>nbkg","NormGen_wt*(signal==0)");
    data_tree->Draw("1>>nsgn","NormGen_wt*(signal==1)");
    cout<<" Pseudodata entries "<<ntot->GetBinContent(1)<<" nsgn "<<nsgn->GetBinContent(1)<<" nbkg "<<nbkg->GetBinContent(1)<<endl;
    if (norm_factor>0){
       RooDataSet dataset_temp ("dataset_tmp","dataset_tmp",RooArgSet(dilep_mass,wgt),RooFit::Import(*data_tree));
       cout<<" before weights"<<endl;
       dataset_temp.Print();
       RooDataSet ds_info("ds_info","ds_info",RooArgSet(dilep_mass,wgt),RooFit::Import(*data_tree),RooFit::WeightVar(wgt));
       cout<<" only gen weight"<<endl;
       ds_info.Print();
       RooRealVar* wgt_branch = (RooRealVar*) dataset_temp.addColumn(wgtFunc);
       RooDataSet dataset("dataset","dataset",(&dataset_temp),*(dataset_temp.get()),0,wgt_branch->GetName()) ;
       cout<<" after weights"<<endl;
       dataset.Print();
       return dataset;
    } else{
       RooDataSet ds("dataset","dataset",RooArgSet(dilep_mass,wgt),RooFit::Import(*data_tree),RooFit::WeightVar(wgt));
       ds.Print();
       return ds;
    }
  } else {
     cout<<" Real data entries "<<data_tree->GetEntries()<<endl;
     return RooDataSet("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));
  } 
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// v12 updates:
//////////////////// +ability to use pseudodata as input and weight it
//////////////////// +use pseudodata or data for combine output instead of toy generated from sideband fit


int ZMuE_fit_v12(TString name="bin1_r2", 
  //  TString data_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root",
    TString data_file="pseudo_data_from_MC_v2_r0.root",
    TString sgn_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root",
    TString sgn_mu_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEUP_bdt_v7_signal_mcRun1*.root", 
    TString sgn_mu_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEDOWN_bdt_v7_signal_mcRun1*.root", 
    TString sgn_ele_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEUP_bdt_v7_signal_mcRun1*.root", 
    TString sgn_ele_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEDOWN_bdt_v7_signal_mcRun1*.root", 
    TString xgbmin="0.3",TString xgbmax="0.7", bool plot_primitives=true,
    bool unblind=false, bool create_dc_input=false, float expected_Nsgn=0,
    TString outvar="mass_ll", bool syst_sgn=false, bool altfit_bkg=true, 
    TString varname="bin", bool pseudodata_input=false, 
    float pseudodata_norm=-1.0, int ntoys=-1, float ntoy_bins=0.5){

   gROOT->SetBatch(true);
   TString dilep_var_name="mass_ll";
   TString dilep_muup_name="mass_ll_Muon_scale_up";
   TString dilep_mudown_name="mass_ll_Muon_scale_down";
   TString dilep_eleup_name="mass_ll_Electron_scale_up";
   TString dilep_eledown_name="mass_ll_Electron_scale_down";
   std::vector<TString> syst_names{"Muon_scale_up", "Muon_scale_down", "Electron_scale_up", "Electron_scale_down"};
   std::vector<TString> syst_files{sgn_mu_up_file, sgn_mu_down_file, sgn_ele_up_file, sgn_ele_down_file};
  
   TString data_combine_file="pseudo_data_from_MC_v2_r0.root";
   bool pseudodata_combine=true;

   TString pseudodata_r0_for_pull="pseudo_data_from_MC_v2_r0.root";
   float signal_toy_r=0.0;  
   bool Save_fittoy_plots=false;
   bool Print_fittoy_pull=true;
 

   bool altfit_bkg_bst3=false;
   bool altfit_bkg_bst4=false;
   bool altfit_bkg_cheb4=false;
   bool altfit_bkg_gamma=false;
   bool altfit_bkg_exp=false; 
 
   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron && mass_ll>70 && mass_ll<110";
   TString tree_name="mytreefit";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;
   double blind_min = 86.;
   double blind_max = 96.;
   int nbin_data = 40;
   bool Verbose=false;
   bool Bkg_only_fit_whole_region=false;


   
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ZMuE_fit: v12 ***Main bkg: Chebychev3 ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;
   //   cout<<" ********************************************************** "<<endl;
   
   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   
   // read trees
   TTree * sgn_tree = get_tree("mytreefit",sgn_file,cuts);
  
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooDataSet sgn_dataset("sgn_dataset","sgn_dataset",RooArgSet(dilep_mass),RooFit::Import(*sgn_tree));
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");

   /////////////////////////////// Signal Fit ////////////////////////////////
   cout<<"\n ************** Signal fit: ************** "<<endl;
   cout<<" Signal MC events: "<< sgn_tree->GetEntries()<<" norm to "<<expected_Nsgn<<endl;
   
   RooWorkspace ws ("ws");
   ws.import(dilep_mass);

   ///// parameters
   RooRealVar sgn_cb_mean("sgn_cb_mean_"+varname,"",91.0e+00, 84.0e+00, 100.0e+00);
   RooRealVar sgn_cb_width("sgn_cb_width_"+varname,"",5., 0.1, 10.);
   RooRealVar sgn_cb_a1("sgn_cb_a1_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n1("sgn_cb_n1_"+varname,"",1.0, 0.00001, 100.0);
   RooRealVar sgn_cb_a2("sgn_cb_a2_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n2("sgn_cb_n2_"+varname,"",1.0, 0.1, 100.0);

   ///// function
   RooDoubleCB sgn_PDF("sgn_PDF","cb",dilep_mass,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2); 

   ///// extend
   RooRealVar n_sgn("sgn_PDF_norm", "",1000,-10,10e10);
   RooAddPdf esgn_PDF("sgn","esgnPDF",RooArgList(sgn_PDF),RooArgList(n_sgn));

   ///// fit
   RooFitResult * sgn_result = esgn_PDF.fitTo(sgn_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range, max_fit_range),RooFit::PrintLevel(-1));
   ///// fix shape
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
      save_plot(sgn_frame,"m(#mu,e)","v12_prmtv_sgn_"+name);
   
   cout<<" Result  whole range    |    85-95 only"<<endl;
   std::pair<double,double> nSgn = yield_calc( n_sgn.getVal(), dilep_mass, &sgn_PDF);
   cout<<" - Expected nSgn "<<nSgn.first<<"   |  "<<nSgn.second<<endl;
   cout<<" ************************************************* "<<endl;
    
   /////// signal systematics
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

   /////// create signal workspace
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
      wspace_sgn->writeToFile("workspace_v12_sgn_"+name+".root");
   }
   
   ///////////////////////////////// Data/bkg fit ////////////////////////////
   cout<<"\n *********************** Data fit ********************** "<<endl;
   TTree * data_tree = get_tree("mytreefit",data_file,cuts);
   
   RooRealVar NormGen_wt("NormGen_wt","gen weight", 0,-10,10);
   RooFormulaVar wgtFunc("wgtFunc","weight formula","NormGen_wt*"+TString(std::to_string(pseudodata_norm)),NormGen_wt);
   RooDataSet dataset = GetDataSet(data_tree, dilep_mass, pseudodata_input, NormGen_wt,wgtFunc,pseudodata_norm);

   cout<<"dataset = "<<dataset.sumEntries()<<endl;
   if (dataset.isWeighted()) cout<<"Is weighted"<<endl;
   else cout<<"Is NOT weighted"<<endl;
   
   ws.import(dataset); 

   dilep_mass.setRange("left",min_fit_range, blind_min); //70,88
   dilep_mass.setRange("right",blind_max, max_fit_range); //94,110
   dilep_mass.setRange("sr",blind_min,blind_max); //88,94
   dilep_mass.setRange("full",min_fit_range, max_fit_range);

   std::vector<RooAbsPdf*> bkg_pdfs;
   std::vector<RooRealVar*> bkg_ampl;
   std::vector<TString> bkg_fnc_names;
   std::vector<TString> bkg_fnc_legs;

   /////// cheb 3rd order (main)
   std::vector<RooRealVar> bkg_cheb3_params = ChebParams(3, varname);
   RooChebychev cheb3_bkgPDF("cheb3_bkgPDF","",dilep_mass,RooArgList(bkg_cheb3_params[0],bkg_cheb3_params[1],bkg_cheb3_params[2]));
   RooRealVar cheb3_n_bkg("cheb3_bkgPDF_"+varname+"_norm","",dataset.numEntries(),0,2*dataset.numEntries());
   bkg_pdfs.push_back(&cheb3_bkgPDF);
   bkg_ampl.push_back(&cheb3_n_bkg);
   bkg_fnc_names.push_back("cheb3_"+name);
   bkg_fnc_legs.push_back("Chebychev 3");

   /////// cheb 4th order
   std::vector<RooRealVar> bkg_cheb4_params = ChebParams(4, varname);
   RooChebychev cheb4_bkgPDF("cheb4_bkgPDF","",dilep_mass,RooArgList(bkg_cheb4_params[0],bkg_cheb4_params[1],bkg_cheb4_params[2],bkg_cheb4_params[3]));
   RooRealVar cheb4_n_bkg("cheb4_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooChebychev cheb4_bkgPDF_out("cheb4_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb4_params[0],bkg_cheb4_params[1],bkg_cheb4_params[2],bkg_cheb4_params[3]));
   if (altfit_bkg && altfit_bkg_cheb4){
      bkg_pdfs.push_back(&cheb4_bkgPDF);
      bkg_ampl.push_back(&cheb4_n_bkg);
      bkg_fnc_names.push_back("cheb4_"+name);
      bkg_fnc_legs.push_back("Chebychev 4");
   }

   /////// 3rd order B. 
   std::vector<RooRealVar> bkg_bst3_params= BstParams(3, varname);
   RooBernsteinFast<3> bst3_bkgPDF("bst3_bkgPDF_"+varname, "", dilep_mass, RooArgList(bkg_bst3_params[0],bkg_bst3_params[1],bkg_bst3_params[2]) );
   bst3_bkgPDF.protectSubRange(true);
   RooRealVar bst3_n_bkg("bst3_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());   
   if (altfit_bkg && altfit_bkg_bst3){
      bkg_pdfs.push_back(&bst3_bkgPDF);
      bkg_ampl.push_back(&bst3_n_bkg);
      bkg_fnc_names.push_back("bst3_"+name);
      bkg_fnc_legs.push_back("Bernstein 3");
    } 

   /////// 4th order B. 
   std::vector<RooRealVar> bkg_bst4_params= BstParams(4, varname);
   RooBernsteinFast<4> bst4_bkgPDF("bst4_bkgPDF", "", dilep_mass, RooArgList(bkg_bst4_params[0],bkg_bst4_params[1],bkg_bst4_params[2],bkg_bst4_params[3]) );
   bst4_bkgPDF.protectSubRange(true);
   RooRealVar bst4_n_bkg("bst4_bkgPDF_norm","",dataset.numEntries(),0,2*dataset.numEntries());   
   if (altfit_bkg && altfit_bkg_bst4){
      bkg_pdfs.push_back(&bst4_bkgPDF);
      bkg_ampl.push_back(&bst4_n_bkg);
      bkg_fnc_names.push_back("bst4_"+name);
      bkg_fnc_legs.push_back("Bernstein 4");
    }

   /////// gamma
   RooRealVar bkg_gamma_a("bkg_gamma_g_"+varname, "",0.7, 0.01, 20.); //5., 0.01, 10.); // 60., 50., 190.
   RooRealVar bkg_gamma_b("bkg_gamma_b_"+varname, "", 26, 0.01, 1500.0); //10, 0.01, 1500.0); // 0.05, 0.01, 0.1
   RooRealVar bkg_gamma_mu("bkg_gamma_mu_"+varname, "", -100, -1000, 70.); //1, -1000, 70.);
   RooGamma gamma_bkgPDF("gamma_bkgPDF","", dilep_mass, bkg_gamma_a,bkg_gamma_b, bkg_gamma_mu);
   RooRealVar gamma_n_bkg("gamma_bkgPDF_norm","",dataset.numEntries()/2,0,2*dataset.numEntries());
   RooGamma gamma_bkgPDF_out("gamma_bkgPDF_"+varname,"", dilep_mass_out, bkg_gamma_a,bkg_gamma_b, bkg_gamma_mu);
    if (altfit_bkg && altfit_bkg_gamma){
       bkg_pdfs.push_back(&gamma_bkgPDF);
       bkg_ampl.push_back(&gamma_n_bkg);
       bkg_fnc_names.push_back("gamma_"+name);
       bkg_fnc_legs.push_back("Gamma");
    }	 

   //////// pol exp
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
   RooAddPdf polexp_bkgPDF("polexp_bkgPDF","",RooArgList(polexp_bkgPDF1,polexp_bkgPDF2),RooArgList(bkg_polexp_c0));
   
   RooRealVar polexp_n_bkg("polexp_bkgPDF_"+varname+"_norm","",dataset.numEntries()/2.,0,1.5*dataset.numEntries());
   
   if (altfit_bkg && altfit_bkg_exp){ 
      bkg_pdfs.push_back(&polexp_bkgPDF);
      bkg_ampl.push_back(&polexp_n_bkg);
      bkg_fnc_names.push_back("expo_"+name);
      bkg_fnc_legs.push_back("#Sigma Expo");
   }

   /////////////////////////////// Fits ///////////////////////////////////
   //////// Bkg only
   FitBkgFunctions(bkg_pdfs,bkg_ampl,dataset,dilep_mass,bkg_fnc_names,bkg_fnc_legs,Bkg_only_fit_whole_region, unblind, nbin_data,pseudodata_input, name);
  // std::vector<double> final_params_cheb3 ={bkg_cheb3_params[0].getVal(),bkg_cheb3_params[1].getVal(),bkg_cheb3_params[2].getVal()};

   ///////// unblinding
   if (unblind)
     FitTotalFunctions(bkg_pdfs, &sgn_PDF, dilep_mass, dataset, bkg_fnc_names, nbin_data, pseudodata_input, name);
   
    
   //////////////////// create bkg and data workspace ////////////////////////
   if (create_dc_input){
     ////// create toy as combine data
     /// bkg toy
     TTree * dc_tree = get_tree("mytreefit",data_combine_file,cuts);
     TH1F* hout = new TH1F("hout","",(max_fit_range-min_fit_range)/0.1,min_fit_range,max_fit_range);
     dc_tree->Draw("mass_ll>>hout","NormGen_wt");
     RooDataHist rd_hout("rd_hout","rd_hout",RooArgSet(dilep_mass_out),hout);
     RooHistPdf pdf_hout("pdf_hout","pdf_hout",RooArgSet(dilep_mass_out),rd_hout); 
     RooDataSet * data_obs = pdf_hout.generate(RooArgSet(dilep_mass_out),dataset.sumEntries()); 
     data_obs->SetName("sim_data_obs_"+varname);
     dilep_mass_out.setBins((max_fit_range-min_fit_range)/0.25);
     RooDataHist hdata_sim("sim_binned_obs_"+varname,"",dilep_mass_out,*data_obs);
     hdata_sim.SetName("sim_binned_obs_"+varname);

     auto dobs_frame = dilep_mass_out.frame();
     data_obs->plotOn(dobs_frame,RooFit::Binning(40),RooFit::MarkerColor(1),RooFit::LineColor(1));
     save_plot(dobs_frame,"m(e,#mu)","v12_data_obs_"+name);

     /// signal toy
     RooDataSet * signal_obs_r1 = sgn_PDF_out.generate(RooArgSet(dilep_mass_out),expected_Nsgn);
     RooDataSet *data_obs_r1 =(RooDataSet*) data_obs->Clone("sim_data_obs_r1_"+varname);
     data_obs_r1->append( *signal_obs_r1);

     ///// output bkg functions
     RooArgList models_out;
     DefaultVals(bkg_cheb3_params,{-1,0.6,-0.2});
     RooChebychev cheb3_bkgPDF_out("cheb3_bkgPDF_"+varname,"",dilep_mass_out,RooArgList(bkg_cheb3_params[0],bkg_cheb3_params[1],bkg_cheb3_params[2]));
     models_out.add(cheb3_bkgPDF_out);

     DefaultVals({bkg_polexp_x0,bkg_polexp_x1,bkg_polexp_c0},{-1.8357e-05,-9.0861e-02,0.5});
     RooExponential polexp_bkgPDF1_out("polexp_bkgPDF1_"+varname,"", dilep_mass_out,bkg_polexp_x0 );
     RooExponential polexp_bkgPDF2_out("polexp_bkgPDF2_"+varname,"", dilep_mass_out,bkg_polexp_x1 );
     RooAddPdf polexp_bkgPDF_out("polexp_bkgPDF_"+varname,"",RooArgList(polexp_bkgPDF1_out,polexp_bkgPDF2_out),RooArgList(bkg_polexp_c0));

     RooWorkspace *wspace_bkg = new RooWorkspace("workspace_background","workspace_background");
     wspace_bkg->import(cheb3_bkgPDF_out);
     wspace_bkg->import(cheb3_n_bkg);

     if (altfit_bkg && altfit_bkg_exp){
        models_out.add(polexp_bkgPDF_out);
        wspace_bkg->import(polexp_bkgPDF_out);
        wspace_bkg->import(polexp_n_bkg);
     }
     RooCategory cat("pdfindex_"+varname, ""); 
     RooMultiPdf multipdf("multipdf_"+varname, "", cat, models_out);
     RooRealVar norm_out("multipdf_"+varname+"_norm","",dataset.numEntries(),0,2*dataset.numEntries());  

     wspace_bkg->import(cat, RooFit::RecycleConflictNodes());
     wspace_bkg->import(multipdf, RooFit::RecycleConflictNodes());
     wspace_bkg->import(norm_out);
     wspace_bkg->import(cheb3_bkgPDF_out);
     wspace_bkg->import(cheb3_n_bkg);   
   
     wspace_bkg->import(*data_obs);
     wspace_bkg->import(*data_obs_r1);
     wspace_bkg->import(hdata_sim);
     wspace_bkg->writeToFile("workspace_v12_bkg_"+name+".root");
   }
    
    
   ///////////////////////////////////  toy tests /////////////////////////////
   if (ntoys>0){
     cout<<"pull test"<<endl;
     TTree * toy_tree = get_tree("mytreefit",pseudodata_r0_for_pull,cuts );
     TH1F* htoy = new TH1F("htoy","",(max_fit_range-min_fit_range)*ntoy_bins,min_fit_range,max_fit_range);
     toy_tree->Draw("mass_ll>>htoy","NormGen_wt*"+TString(to_string(pseudodata_norm)));
     TCanvas * c1_toy = new TCanvas("c1_toy","",700,700);
     htoy->Draw();
     c1_toy->SaveAs("hgentoy_"+name+".png");
     htoy->Scale(1./htoy->Integral());
     
     RooDataHist rd_htoy("rd_htoy","rd_htoy",RooArgSet(dilep_mass_out),htoy);
     RooHistPdf pdf_htoy("pdf_htoy","pdf_htoy",RooArgSet(dilep_mass_out),rd_htoy);
     TH1F * hpull_cheb3 = new TH1F("hpull_cheb3","",100,-5,5);
     TH1F * hpull_exp = new TH1F("hpull_exp","",100,-5,5);
     for(int itoy=0; itoy<ntoys; itoy++){
       if (itoy%(int(0.1*ntoys))==0) cout<<" generate toy "<<itoy<<" / "<<ntoys<<endl;
       TString str_toy =TString(to_string(itoy));
       RooDataSet * dataset_toy = pdf_htoy.generate(RooArgSet(dilep_mass_out),dataset.sumEntries());
       RooDataSet * signal_dataset_toy = sgn_PDF_out.generate(RooArgSet(dilep_mass_out),expected_Nsgn*signal_toy_r);
       if (signal_toy_r>0)
          dataset_toy->append( *signal_dataset_toy);
       DefaultVals(bkg_cheb3_params,{-1,0.6,-0.2});
       RooChebychev cheb3_bkgPDF_toy("cheb3_bkgPDF_toy_"+str_toy,"",dilep_mass_out,RooArgList(bkg_cheb3_params[0],bkg_cheb3_params[1],bkg_cheb3_params[2]));
       RooRealVar nBkg_cheb3_toy("nbkg_cheb3_toy","",dataset.sumEntries(),0,2*dataset.sumEntries()); 
       RooRealVar nSgn_cheb3_toy("nsgn_cheb3_toy","",0,-expected_Nsgn*100,expected_Nsgn*100);
       RooAddPdf total_pdf_cheb3_toy("total_pdf_cheb3_"+str_toy,"", RooArgList(cheb3_bkgPDF_toy, sgn_PDF_out),  RooArgList(nBkg_cheb3_toy,nSgn_cheb3_toy));
       RooFitResult* toy_cheb3_result = total_pdf_cheb3_toy.fitTo(*dataset_toy,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range , max_fit_range),RooFit::PrintLevel(-1));
       hpull_cheb3->Fill((nSgn_cheb3_toy.getVal()-expected_Nsgn*signal_toy_r)/nSgn_cheb3_toy.getError());

       DefaultVals({bkg_polexp_x0,bkg_polexp_x1,bkg_polexp_c0},{-1.8357e-05,-9.0861e-02,0.5});
       RooExponential polexp_bkgPDF1_toy("polexp_bkgPDF1_toy_"+str_toy,"", dilep_mass_out,bkg_polexp_x0 );
       RooExponential polexp_bkgPDF2_toy("polexp_bkgPDF2_toy_"+str_toy,"", dilep_mass_out,bkg_polexp_x1 );
       RooAddPdf polexp_bkgPDF_toy("polexp_bkgPDF_toy_"+str_toy,"",RooArgList(polexp_bkgPDF1_toy,polexp_bkgPDF2_toy),RooArgList(bkg_polexp_c0));
       RooRealVar nBkg_exp_toy("nbkg_exp_toy","",dataset.sumEntries(),0,2*dataset.sumEntries());
       RooRealVar nSgn_exp_toy("nsgn_exp_toy","",0,-expected_Nsgn*100,expected_Nsgn*100);
       RooAddPdf total_pdf_exp_toy("total_pdf_exp_"+str_toy,"", RooArgList(polexp_bkgPDF_toy, sgn_PDF_out),  RooArgList(nBkg_exp_toy,nSgn_exp_toy));
       RooFitResult* toy_exp_result = total_pdf_exp_toy.fitTo(*dataset_toy,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range , max_fit_range),RooFit::PrintLevel(-1));
       if (altfit_bkg && altfit_bkg_exp)
          hpull_exp->Fill((nSgn_exp_toy.getVal()-expected_Nsgn*signal_toy_r)/nSgn_exp_toy.getError());
       
       
       if (Print_fittoy_pull){ 
          cout<<" Cheb3: #toy "<<itoy<<" pull: "<<(nSgn_cheb3_toy.getVal()-expected_Nsgn*signal_toy_r)/nSgn_cheb3_toy.getError()<<" Fit: nSgn "<<nSgn_cheb3_toy.getVal()<<" +/- "<<nSgn_cheb3_toy.getError()<<" nBkg "<<nBkg_cheb3_toy.getVal()<<" +/- "<<nBkg_cheb3_toy.getError()<<" entries "<<dataset.sumEntries()<<"; Gen Sgn "<<expected_Nsgn*signal_toy_r<<endl;
          if (altfit_bkg && altfit_bkg_exp)
             cout<<" Exp: #toy "<<itoy<<" pull: "<<(nSgn_exp_toy.getVal()-expected_Nsgn*signal_toy_r)/nSgn_exp_toy.getError()<<" Fit: nSgn "<<nSgn_exp_toy.getVal()<<" +/- "<<nSgn_exp_toy.getError()<<" nBkg "<<nBkg_exp_toy.getVal()<<" +/- "<<nBkg_exp_toy.getError()<<" entries "<<dataset.sumEntries()<<"; Gen Sgn "<<expected_Nsgn*signal_toy_r<<endl;
       }
       if (Save_fittoy_plots){
         auto tmp_toy_frame = dilep_mass_out.frame();
         dataset_toy->plotOn(tmp_toy_frame,RooFit::Binning(nbin_data));
         total_pdf_cheb3_toy.plotOn(tmp_toy_frame,RooFit::LineColor(kBlack),RooFit::Range("full"),RooFit::NormRange("full"),RooFit::Name("total_pdf_"+str_toy));
         cheb3_bkgPDF_toy.plotOn(tmp_toy_frame,RooFit::LineColor(kRed),RooFit::Normalization(nBkg_cheb3_toy.getVal(), RooAbsReal::NumEvent),RooFit::Name("background_cheb3_pdf_"+str_toy),RooFit::LineStyle(7));
         sgn_PDF_out.plotOn(tmp_toy_frame,RooFit::LineColor(kBlue),RooFit::Normalization(nSgn_cheb3_toy.getVal(), RooAbsReal::NumEvent),RooFit::Name("signal_pdf_"+str_toy),RooFit::LineStyle(7));
         save_plot(tmp_toy_frame,"m(#mu,e)","v12_fit_toy_cheb3_"+str_toy);   
         auto tmp_toy_frame2 = dilep_mass_out.frame();
         dataset_toy->plotOn(tmp_toy_frame2,RooFit::Binning(nbin_data));
         total_pdf_exp_toy.plotOn(tmp_toy_frame2,RooFit::LineColor(kBlack),RooFit::Range("full"),RooFit::NormRange("full"),RooFit::Name("total_pdf_"+str_toy));
         polexp_bkgPDF_toy.plotOn(tmp_toy_frame2,RooFit::LineColor(kRed),RooFit::Normalization(nBkg_exp_toy.getVal(), RooAbsReal::NumEvent),RooFit::Name("background_exp_pdf_"+str_toy),RooFit::LineStyle(7));
         sgn_PDF_out.plotOn(tmp_toy_frame2,RooFit::LineColor(kBlue),RooFit::Normalization(nSgn_exp_toy.getVal(), RooAbsReal::NumEvent),RooFit::Name("signal_pdf_"+str_toy),RooFit::LineStyle(7));
         if (altfit_bkg && altfit_bkg_exp)
            save_plot(tmp_toy_frame2,"m(#mu,e)","v12_fit_toy_exp_"+str_toy); 
       }
       
     }  
     TCanvas * cpull = new TCanvas("cpull","",700,700);
     hpull_cheb3->Draw();
     hpull_exp->SetLineColor(2);
     if (altfit_bkg && altfit_bkg_exp)
        hpull_exp->Draw("sames");
     cpull->SaveAs("pull_"+name+"_toyR"+TString(to_string(signal_toy_r))+".png");    
   }


 return 0;
} 
