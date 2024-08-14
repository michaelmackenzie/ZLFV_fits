#ifndef FIT_HELPER_CORE_H
#define FIT_HELPER_CORE_H

#include <fstream>
#include <random>
#include "fit_helper.h"

/////////////////////////// Core helper functions ////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// used from v15
//////////////////// 
//////////// update 1/4/24: GetDataSet -> reads also th1 and generate toy

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Input/generic //////////////////////////////
///////////////////////////////////////////////////////////////////////////////



///////// get dataset from data/MC and rescales the latter with weight
RooDataSet GetDataSet(TString data_file, TString cuts, RooRealVar &dilep_mass, bool Is_pseudodata, RooRealVar &wgt, RooFormulaVar &wgtFunc, float norm_factor, bool Is_histo){
  
  if (Is_pseudodata){
    TTree * data_tree = get_tree("mytreefit",data_file,cuts);
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
  } else if (Is_histo){
    TFile * fhisto = new TFile(data_file,"READ");
    TH1F * histo = (TH1F*) fhisto->Get("hbkg_safe");
    RooDataHist data_histo("data_histo","data_histo",RooArgSet(dilep_mass),histo);
    RooHistPdf pdf_histo("pdf_histo","",RooArgSet(dilep_mass),data_histo);
    RooDataSet * sim_dataset = pdf_histo.generate(RooArgSet(dilep_mass),histo->Integral());
    sim_dataset->SetName("dataset");
    return *sim_dataset;
  } else {
     TTree * data_tree = get_tree("mytreefit",data_file,cuts);
     cout<<" Real data entries "<<data_tree->GetEntries()<<endl;
     return RooDataSet("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));
  } 
}


///////// get dataset from data only - overloaded simpler function
RooDataSet GetDataSet(TTree *data_tree, RooRealVar &dilep_mass){
  cout<<" Real data entries "<<data_tree->GetEntries()<<endl;
  return RooDataSet("dataset","dataset",RooArgSet(dilep_mass),RooFit::Import(*data_tree));
   
}


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

//////// generic plot function 
void PlotFunctions(std::vector<RooAbsPdf*> pdfs, RooPlot * xframe, RooRealVar dilep_mass, std::vector<TString> names, std::vector<TString> legs, TLegend *leg, RooDataSet dataset, TString name, int nbin_data, bool unblind){
   TString norm_range="left,right";
   if (unblind)
      norm_range="full";
   for (int i=0; i<pdfs.size(); i++){
     pdfs[i]->plotOn(xframe,RooFit::LineColor(i+1),RooFit::Range("full"),RooFit::NormRange(norm_range),RooFit::Name(names[i]));
     save_pull(xframe, dilep_mass, "m(e,#mu)","v15_"+names[i]);
     leg->AddEntry(xframe->findObject(names[i]),legs[i]);
   }
   if (unblind)
     dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::Name("data"));
   else
     dataset.plotOn(xframe,RooFit::Binning(nbin_data),RooFit::CutRange("left,right"),RooFit::Name("data"));
   leg->AddEntry(xframe->findObject("data"),"Data");
   save_plot(xframe,"m(e,#mu)","v15_sr_total_"+name,leg);
}



///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Signal part ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



///////// systematic functions parameters : returns maximum deviation of mean/width
std::pair<double,double> SignalSystematicsMaxMeanWidth(TString syst_file, TString cuts, TString syst_name, float min_fit_range, float max_fit_range, int nbin_data,TString outname){
  
    TTree * syst_tree = get_tree("mytreefit",syst_file,cuts); 
    RooRealVar dilep_mass_syst("mass_ll_"+syst_name,"m(e,#mu)", (min_fit_range-max_fit_range)/2., min_fit_range , max_fit_range, "GeV/c^{2}");
    RooDataSet syst_dataset("syst_"+syst_name+"_dataset","",RooArgSet(dilep_mass_syst),RooFit::Import(*syst_tree));
    RooRealVar syst_cb_mean("sgn_cb_mean_"+syst_name,"",91.0e+00, 84.0e+00, 100.0e+00);
    RooRealVar syst_cb_width("syst_cb_width_"+syst_name,"",2.5, 0.1, 10.);
    RooRealVar syst_cb_a1("syst_cb_a1_"+syst_name,"",1.0, 0.00001, 10.0);
    RooRealVar syst_cb_n1("syst_cb_n1_"+syst_name,"",3.0, 0.00001, 10.0);
    RooRealVar syst_cb_a2("syst_cb_a2_"+syst_name,"",2.5, 0.00001, 10.0);
    RooRealVar syst_cb_n2("syst_cb_n2_"+syst_name,"",1.3, 0.00001, 10.0);

    RooDoubleCB syst_sgn_PDF("sgn_PDF_"+syst_name+"_syst","cb",dilep_mass_syst,syst_cb_mean,syst_cb_width,syst_cb_a1,syst_cb_n1,syst_cb_a2,syst_cb_n2);
    RooRealVar n_syst("n_syst"+syst_name, "",1000,0,10000000000);
    RooAddPdf syst_esgn_PDF("esgn_PDF_"+syst_name,"",RooArgList(syst_sgn_PDF),RooArgList(n_syst));
    // fit result
    RooFitResult * syst_result = syst_esgn_PDF.fitTo(syst_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range , max_fit_range),RooFit::PrintLevel(-1));

    double mean=syst_cb_mean.getVal(), width=syst_cb_width.getVal();
    auto syst_frame = dilep_mass_syst.frame();
    print_details (syst_result);
    int n_param_syst = syst_result->floatParsFinal().getSize();
    dilep_mass_syst.setBins(nbin_data);
    float chi2_syst = get_chi_squared(dilep_mass_syst, &syst_esgn_PDF, syst_dataset, true, nbin_data, n_param_syst); 
    syst_dataset.plotOn(syst_frame,RooFit::Binning(nbin_data));
    syst_esgn_PDF.plotOn(syst_frame,RooFit::LineColor(kBlue));
    save_plot(syst_frame,"m(#mu,e)","v15_prmtv_syst_"+syst_name+"_"+outname);
    delete syst_result;
    return std::make_pair(mean,width);
}





///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// BKG PDFs ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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


/////////// chebychev
RooChebychev * CreateChebychev( TString name, int order, RooRealVar& dilep_mass, std::vector<float> parameters={}){
   float def=0,min=-1.0,max=1.0;
   RooArgList param_list;
   std::vector<RooRealVar*> bkg_params;
   for (int i=0; i<order; i++){
     if (i==0)
        def=-1,min=-5.0,max=5.0;
     else if (i==1) def=0.6;
     else if (i==2) def=-0.2;
     if (parameters.size()>0) def= parameters[i];
     bkg_params.push_back(new RooRealVar(name+"_"+TString(to_string(i)), "", def,min,max));
     param_list.add(*bkg_params[i]);
   }
  return new  RooChebychev(name,"",dilep_mass,param_list);
}

//////////// sum of expo
RooAddPdf * CreateSumExpo(TString name, int order, RooRealVar& dilep_mass, bool recursive_coef=true, std::vector<float> parameters={})
{
   std::vector<RooRealVar*> bkg_params;  
   std::vector<RooExponential*> exps;
   std::vector<RooRealVar*> coefs;
   RooArgList exp_list;
   RooArgList coef_list;

   float def=-0.1,min=-10,max=10;
   for( int i=0; i<order; i++){
     if (order ==1 ) def=-0.07;
     if (order ==2 && i==0 ) def=0.05,max=10;
     if (order ==2 && i==1 ) def=-0.1,max=10;
     if (order >2 && i==0 ) def=0.05,max=10;
     if (order >2 && i==1 ) def=-0.1,max=10;
     if (order >2 && i>1 ) def=-1*i,max=0;
     if (parameters.size()>0 && recursive_coef) def= parameters[i+order];
     if (parameters.size()>0 && !recursive_coef) def= parameters[i+order+1];
     RooRealVar * bkg_sumexp_x = new RooRealVar(name+"_x"+TString(to_string(i)), "",def,min,max);
     bkg_params.push_back(bkg_sumexp_x);
     RooExponential * sumexp_bkgPDF = new RooExponential(name+"_bkgPDF"+TString(to_string(i)),"", dilep_mass, *bkg_sumexp_x );
     exps.push_back(sumexp_bkgPDF);
     exp_list.add(*exps[i]);
     if (i>0 && recursive_coef){
        float c_def = 1./order;
        if (parameters.size()>0) c_def= parameters[i-1];
        RooRealVar * bkg_sumexp_c = new RooRealVar(name+"_c"+TString(to_string(i)), "",c_def, 0, 1.);
        coefs.push_back(bkg_sumexp_c);
        coef_list.add(*coefs[i-1]);
     } 
    if (!recursive_coef){
        float c_def = 1./order;
        if (parameters.size()>0) c_def= parameters[i];
        RooRealVar * bkg_sumexp_c = new RooRealVar(name+"_c"+TString(to_string(i)), "",c_def, 0., 100.);
        coefs.push_back(bkg_sumexp_c);
        coef_list.add(*coefs[i]);
     }
     
 
   }

   return new RooAddPdf(name,"",exp_list,coef_list,recursive_coef);
}


//////////// power law
RooAddPdf * CreateSumPower(TString name, int order, RooRealVar& dilep_mass,  bool recursive_coef=true,  std::vector<float> parameters={})
{
   std::vector<RooRealVar*> bkg_params;
   std::vector<RooPower*> powers;
   std::vector<RooRealVar*> coefs;
   RooArgList power_list;
   RooArgList coef_list;
   RooArgSet coef_set;

   float def=0,min=-200,max=200;
   for( int i=0; i<order; i++){
     if (order==1) def=-10;
     if (order==2 && i==0) def=10;
     if (order==2 && i==1) def=-8;
     if (order>2 && i==0) def=-10,min=-100,max=0;
     if (order>2 && i==1) def=-8,min=-100,max=0;
     if (order>2 && i==2) def=20,min=-100,max=100;

     if (parameters.size()>0) def= parameters[i]; 
     RooRealVar * bkg_sumplaw_a = new RooRealVar(name+"_a"+TString(to_string(i)), "",def,min,max);
     bkg_params.push_back(bkg_sumplaw_a);
     
     RooPower * sumplaw_bkgPDF = new RooPower(name+"_bkgPDF"+TString(to_string(i)), "", dilep_mass, *bkg_sumplaw_a);
     powers.push_back(sumplaw_bkgPDF);
     power_list.add(*powers[i]);
    
     if ( i<order-1 && recursive_coef){
        float c_def = 1./order;
        if (parameters.size()>0) c_def= parameters[i+order];
        RooRealVar * bkg_sumplaw_c = new RooRealVar(name+"_c"+TString(to_string(i)), "",c_def, 0, 1.); 
        coefs.push_back(bkg_sumplaw_c);
        coef_list.add(*coefs[i]);
     }
     
     if ( !recursive_coef){
        float c_def = 1./order;
        if (parameters.size()>0) c_def= parameters[i+order];
        RooRealVar * bkg_sumplaw_c = new RooRealVar(name+"_c"+TString(to_string(i)), "",c_def, 0, 1.);
        coefs.push_back(bkg_sumplaw_c);
        coef_list.add(*coefs[i]);
     }
 
   }
    
   return new RooAddPdf(name,"",power_list,coef_list);
}






///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Fits ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////// Fit bkg only functions sideband or total fit and blind/unblind plot
std::vector<std::vector<float>> FitBkgFunctions(std::vector<RooAbsPdf*> pdfs, std::vector<RooRealVar*> ampls, RooDataSet &dataset, RooRealVar &dilep_mass, vector<TString> names, std::vector<TString> legs, bool Bkg_only_fit_whole_region, bool unblind, int nbin_data, int nbin_blind, bool is_pseudodata , TString extra_name, bool Print_details=true){


 std::vector<std::vector<float>> output;
 TString fit_range="left,right";
  if (Bkg_only_fit_whole_region) 
     fit_range="full";
  TString plot_range="left,right";
  if (unblind)
     plot_range="full";
  
  TLegend * leg = new TLegend(0.4,0.9-names.size()*0.1,0.9,0.9);
  auto plot_frame =  dilep_mass.frame();
  TPaveText * pt = new TPaveText(0.4,0.9-names.size()*0.1,0.9,0.9,"tlNDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetTextAlign(12);
  pt->AddText("#chi^{2}:\n");

  dataset.plotOn(plot_frame,RooFit::Binning(nbin_data),RooFit::Name("data"),RooFit::AsymptoticError(is_pseudodata),RooFit::LineColor(0),RooFit::MarkerColor(0));

  cout<<"\n *********** Bkg-only fit in "+fit_range<<" range *********"<<endl;  
  for (int i=0; i<pdfs.size(); i++){ 
    std::vector<float> temp;
    cout<<" Fit: "+names[i]<<endl;
    RooAddPdf epdf("epdf_"+names[i],"", RooArgList(*pdfs[i]),  RooArgList(*ampls[i]));
   
    RooFitResult * fit_result = epdf.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range(fit_range),RooFit::AsymptoticError(is_pseudodata));
    if (Print_details) print_details (fit_result);
    int n_param = fit_result->floatParsFinal().getSize();
    dilep_mass.setBins(nbin_data);
    
    float chi2 = get_chi_squared(dilep_mass, pdfs[i], dataset, unblind, nbin_data-nbin_blind, n_param, false,is_pseudodata);
    float pvalue = ROOT::Math::chisquared_cdf_c(chi2,nbin_data-nbin_blind-n_param);
    if (chi2<0) chi2=1000000000;
    temp.push_back(chi2);
//    temp.push_back(fit_result->minNll());
    temp.push_back(nbin_data-nbin_blind-n_param);
    
    for (unsigned int ik=0; ik<fit_result->floatParsFinal().getSize(); ik++)
      temp.push_back(static_cast<RooAbsReal*>(fit_result->floatParsFinal().at(ik))->getVal());      
    
    output.push_back(temp);
    cout<<">>>>> p-value "<<pvalue<<endl;
    pdfs[i]->plotOn(plot_frame,RooFit::LineColor(i+1),RooFit::Range("full"),RooFit::Name(names[i]));
    leg->AddEntry(plot_frame->findObject(names[i]), Form(legs[i]+" (p: %1.2lf, #chi^2: %2.4lf)",pvalue, chi2/(nbin_data-nbin_blind-n_param) ));
    pt->AddText(Form(legs[i]+" %3.4lf",pvalue));
    if (Print_details){
      cout<<" ******** Bkg-only result ********* "<<endl;
      cout<<"whole range    |    85-95 only"<<endl;
   }
    std::pair<double,double> nBkg = yield_calc( ampls[i]->getVal(), dilep_mass, pdfs[i]);
   if (Print_details)
      cout<<" - nBkg("+names[i]+") "<<nBkg.first<<"  |  "<<nBkg.second<<endl;    
  }
  
  dataset.plotOn(plot_frame,RooFit::Binning(nbin_data),RooFit::CutRange(plot_range),RooFit::Name("data"),RooFit::SumW2Error(is_pseudodata) );
  save_plot(plot_frame,"m(e,#mu)","v15_bkgfit_"+extra_name,leg);
  cout<<" ***************************************"<<endl;
  return output;
}

/////////// adds bkg and signal functions and fit them unblind ///////////////
std::vector<std::pair<float,int>> FitTotalFunctions(std::vector<RooAbsPdf*> bkg_pdfs, RooAbsPdf* signal_pdf, RooRealVar dilep_mass, RooDataSet &dataset, vector<TString> names, int nbin_data, bool is_pseudodata, TString extra_name){


  cout<<"\n ************ Unblind fit ************** "<<endl;
  std::vector<std::pair<float,int>> output;
  for (int i=0; i<bkg_pdfs.size(); i++){
    cout<<" PDF: "<<names[i]<<endl;
    RooRealVar nBkg_tmp("nBkg_"+names[i],"", dataset.sumEntries(),0.5*dataset.numEntries() , 2*dataset.sumEntries());
    RooRealVar nSgn_tmp("nSgn_"+names[i],"", 0, -dataset.sumEntries()/10 , dataset.sumEntries()/10);
    RooAddPdf total_pdf("total_pdf_"+names[i],"", RooArgList(*bkg_pdfs[i], *signal_pdf),  RooArgList(nBkg_tmp,nSgn_tmp));
    RooFitResult * fit_result = total_pdf.fitTo(dataset,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("full"),RooFit::AsymptoticError(is_pseudodata));
    //,RooFit::SumW2Error(true));
    //,RooFit::AsymptoticError(is_pseudodata));
    auto tmpframe = dilep_mass.frame();
    print_details (fit_result);
    int n_param = fit_result->floatParsFinal().getSize();
    dilep_mass.setBins(nbin_data);
    float chi2 = get_chi_squared(dilep_mass, &total_pdf, dataset, true, nbin_data, n_param,false,is_pseudodata);
    float pvalue = ROOT::Math::chisquared_cdf_c(chi2,nbin_data-n_param);
    output.push_back(std::make_pair(chi2,nbin_data-n_param));
    cout<<">>>>> p-value "<<pvalue<<endl;
    dataset.plotOn(tmpframe,RooFit::Binning(nbin_data),RooFit::Name("data"),RooFit::AsymptoticError(is_pseudodata));
    bkg_pdfs[i]->plotOn(tmpframe,RooFit::LineColor(kRed),RooFit::Normalization(nBkg_tmp.getVal(), RooAbsReal::NumEvent),RooFit::Name("background_pdf_"+names[i]),RooFit::LineStyle(7));
    signal_pdf->plotOn(tmpframe,RooFit::LineColor(kBlue),RooFit::Normalization(nSgn_tmp.getVal(), RooAbsReal::NumEvent),RooFit::Name("signal_pdf_"+names[i]),RooFit::LineStyle(7));
    total_pdf.plotOn(tmpframe,RooFit::LineColor(kBlack),RooFit::NormRange("full"), RooFit::Name("total_pdf_"+names[i]));
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
    
    pt->AddText("p-value: "+TString(std::to_string( pvalue )));
    pt->AddText(Form("#chi^{2}: %1.2f",chi2/(nbin_data-n_param)));
    save_plot(tmpframe,"m(e,#mu)","v15_totalfit_"+names[i],leg,pt);
    delete leg;
  }
  cout<<" ****************************** "<<endl;
  return output;
}


std::pair<int,float> Ftest(std::vector<RooAbsPdf*> pdfs, std::vector<RooRealVar*> ampls, RooDataSet &dataset, RooRealVar &dilep_mass, vector<int> orders, vector<TString> names, std::vector<TString> legs, int nbin_data, int nbin_blind, bool unblind_data_or_ispseudo, bool is_pseudodata , TString extra_name, float ftest_step, float min_pvalue=-1, int print_level=0 ){
  
  bool Print_details = 0;
  if (print_level) Print_details=true;
  std::vector<std::vector<float>> chi2_dof = FitBkgFunctions(pdfs, ampls, dataset, dilep_mass, names, legs, unblind_data_or_ispseudo, unblind_data_or_ispseudo, nbin_data,nbin_blind, is_pseudodata , "ftest_"+extra_name,Print_details);

  std::vector<float> pvalues; 
  std::vector<float> ftests;
  int temp_max_order=0;
  float temp_max_pval=0;
  for (int i =0; i<chi2_dof.size(); i++){
    pvalues.push_back(ROOT::Math::chisquared_cdf_c(chi2_dof[i][0],chi2_dof[i][1]) );
    if (ROOT::Math::chisquared_cdf_c(chi2_dof[i][0],chi2_dof[i][1])>temp_max_pval){
      temp_max_order=orders[i];
      temp_max_pval=ROOT::Math::chisquared_cdf_c(chi2_dof[i][0],chi2_dof[i][1]);
    }  
    if (i>0)
       ftests.push_back( ROOT::Math::chisquared_cdf_c(chi2_dof[i-1][0] - chi2_dof[i][0], chi2_dof[i-1][1]-chi2_dof[i][1]) );
  }  

  vector<int> orders_pass_cut;
  vector<float> chis_pass_cut;
  vector<int> dofs_pass_cut;  
  for (int i =0; i<chi2_dof.size(); i++){
    if (pvalues[i]<min_pvalue) continue;
    orders_pass_cut.push_back(orders[i]);
    chis_pass_cut.push_back(chi2_dof[i][0]);
    dofs_pass_cut.push_back(chi2_dof[i][1]);
  }

  vector<float> ftests_pass_cut;
  for (int i =1; i<orders_pass_cut.size(); i++)
     ftests_pass_cut.push_back( ROOT::Math::chisquared_cdf_c( chis_pass_cut[i-1] - chis_pass_cut[i], dofs_pass_cut[i-1]-dofs_pass_cut[i]) );
     
  int min_order=temp_max_order,best_ftest=0;
  if (chis_pass_cut.size()>0)
     min_order=orders_pass_cut[0];

  for (int i =0; i<ftests_pass_cut.size(); i++){
    if (ftests_pass_cut[i]>ftest_step || orders_pass_cut[i]+1<orders_pass_cut[i+1]) break;
    min_order=orders_pass_cut[i+1];
    best_ftest=ftests_pass_cut[i]; 
  }
  
  TString summary="\n Ftest result";
  for (int i =0; i<pvalues.size(); i++)
    summary+=Form(" order_%1d_chi2_%2.4lf_pvalue_%2.4lf",orders[i],chi2_dof[i][0],pvalues[i]);

  summary+="\n";

  for (int i=0; i<ftests.size(); i++)
    summary+=Form(" order_%1d_to_%1d_ftest_%2.4lf",orders[i],orders[i+1],ftests[i]);

 summary+="\n";

  for (int i=0; i<ftests_pass_cut.size(); i++)
    summary+=Form(" order_%1d_to_%1d_ftest_%2.4lf",orders_pass_cut[i],orders_pass_cut[i+1],ftests_pass_cut[i]);

  cout<<summary<<"; BEST: "<<min_order<<" with "<<best_ftest<<endl;
  
  return std::pair<int,float> (min_order,best_ftest);  
}


///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Output/tests //////////////////////////////
///////////////////////////////////////////////////////////////////////////////


///////////// create data obs from toys 
RooDataSet * GetDataObs(RooRealVar & dilep_mass_out, bool pseudodata_fit_combine, bool data_blinded_combine, bool histo_template_combine, bool pseudodata_template_combine, TString data_combine_file, TString cuts, TString data_histo_file, TString histo_name_combine, float min_fit_range, float max_fit_range, float blind_min, float blind_max,float sumEntries, TString ver, TString name){


  RooDataSet * data_obs=NULL;

  if ( pseudodata_fit_combine || data_blinded_combine ){
     TTree * dc_tree = get_tree("mytreefit",data_combine_file,cuts);
     TH1F* hmc_psd = new TH1F("hmc_psd","",(max_fit_range-min_fit_range),min_fit_range,max_fit_range);
     if (pseudodata_fit_combine)
        dc_tree->Draw("mass_ll>>hmc_psd","NormGen_wt");
     else
        dc_tree->Draw("mass_ll>>hmc_psd","mass_ll<"+TString(std::to_string(blind_min))+" || mass_ll>"+TString(std::to_string(blind_max)));
     RooDataHist rdmc_psd("rdmc_psd","rdmc_psd",RooArgSet(dilep_mass_out),hmc_psd);
     /// fit to get the function params
     std::vector<RooRealVar> psd_cheb3_params= ChebParams(3, "test_psd");
     RooChebychev cheb3_psdPDF("cheb3_psdPDF","",dilep_mass_out,RooArgList(psd_cheb3_params[0],psd_cheb3_params[1],psd_cheb3_params[2]));
     RooRealVar ampl_psd = RooRealVar("ampl_psd","",1000,-100,100000000);
     RooAddPdf epsdPDF("epsdPDF","", RooArgList(cheb3_psdPDF), RooArgList(ampl_psd));
     RooFitResult * epsd_fit_result = NULL;
     if (data_blinded_combine){
        dilep_mass_out.setRange("left",min_fit_range, blind_min);
        dilep_mass_out.setRange("right",blind_max, max_fit_range);
        epsd_fit_result = epsdPDF.fitTo(rdmc_psd,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range("left,right"));
     } else
        epsd_fit_result = epsdPDF.fitTo(rdmc_psd,RooFit::Extended(1),RooFit::Save(),RooFit::PrintLevel(-1));

    /// plot of fit for debug
    auto psd_frame = dilep_mass_out.frame();
    rdmc_psd.plotOn(psd_frame,RooFit::Binning((max_fit_range-min_fit_range)),RooFit::MarkerColor(1),RooFit::LineColor(1));

    if (data_blinded_combine)
       epsdPDF.plotOn(psd_frame,RooFit::MarkerColor(2),RooFit::LineColor(2),RooFit::NormRange("left,right"));
    else
       epsdPDF.plotOn(psd_frame,RooFit::MarkerColor(2),RooFit::LineColor(2),RooFit::Normalization(ampl_psd.getVal(), RooAbsReal::NumEvent));

    save_plot(psd_frame,"m(e,#mu)",ver+"_fit_for_data_obs_"+name);

    for (int i=0; i<3; i++)
      psd_cheb3_params[i].setConstant(true); //set constant to generate toy 

    data_obs = epsdPDF.generate(RooArgSet(dilep_mass_out),sumEntries); //actual generation

  } else if (histo_template_combine){
    TFile * ftemplate = new TFile(data_histo_file,"READ");
    TH1F * htemplate = (TH1F*) ftemplate->Get(histo_name_combine);
    RooDataHist rd_htemplate("rd_htemplate","rd_htemplate",RooArgSet(dilep_mass_out),htemplate);
    RooHistPdf pdf_htemplate("pdf_htemplate","pdf_htemplate",RooArgSet(dilep_mass_out),rd_htemplate);
    data_obs = pdf_htemplate.generate(RooArgSet(dilep_mass_out),sumEntries);

  } else if ( pseudodata_template_combine ){
     TTree * dc_tree = get_tree("mytreefit",data_combine_file,cuts);
     TH1F* hout = new TH1F("hout","",(max_fit_range-min_fit_range),min_fit_range,max_fit_range);
     dc_tree->Draw("mass_ll>>hout","NormGen_wt");
     RooDataHist rd_hout("rd_hout","rd_hout",RooArgSet(dilep_mass_out),hout);
     RooHistPdf pdf_hout("pdf_hout","pdf_hout",RooArgSet(dilep_mass_out),rd_hout);
     data_obs = pdf_hout.generate(RooArgSet(dilep_mass_out),sumEntries);
  }


  return data_obs;

}


TH1F* GetHistoTemplate(TString data_combine_file, TString cuts, float min_fit_range, float max_fit_range, TString histo_name_combine, TString data_histo_file, bool read_th1, TString name="" )
{
   TH1F* hmc_fnc = new TH1F("hmc_fnc_"+name,"",(max_fit_range-min_fit_range),min_fit_range,max_fit_range);

   if (read_th1){
     TFile * ftemplate = new TFile(data_histo_file,"READ");
     hmc_fnc = (TH1F*) ftemplate->Get(histo_name_combine);    
   } else{   
     TTree * dc_tree = get_tree("mytreefit",data_combine_file,cuts);
     dc_tree->Draw("mass_ll>>hmc_fnc_"+name,"NormGen_wt");

   }
   hmc_fnc->Scale(1./hmc_fnc->Integral());
   
    return hmc_fnc;
}



////// generate toys
std::vector<RooDataSet *> GenerateTemplateToys( RooAbsPdf * bkg_pdf, RooAbsPdf * sgn_pdf, RooRealVar & dilep_mass_out, int ntoys, float bkgEntries, float sgnEntries){

  vector<RooDataSet *> datasets;
  for ( int itoy=0; itoy<ntoys; itoy++){
      if (itoy%(int(0.1*ntoys))==0) 
         cout<<" generate toy "<<itoy<<" / "<<ntoys<<endl;
     RooDataSet * bkg_toy = bkg_pdf->generate(RooArgSet(dilep_mass_out), bkgEntries);
     RooDataSet * sgn_toy = sgn_pdf->generate(RooArgSet(dilep_mass_out),sgnEntries);
     if (sgnEntries>0)
        bkg_toy->append( *sgn_toy);
     datasets.push_back(bkg_toy);
  }
  return datasets;  
}



////// test toys
TH1F *  ChebychevPullFromToys( int order, RooAbsPdf * sgn_pdf, RooRealVar & dilep_mass_out, std::vector<RooDataSet *> toys, float nExpected, float sumEntries, TString name,bool PrintPulls=false){
  float range_nSgn = nExpected*100;
  if (nExpected==0) 
     range_nSgn =100;
  TH1F * hpull_cheb = new TH1F("hpull_cheb_"+name,"",100,-5,5);
  for (int itoy=0; itoy<toys.size(); itoy++){
     if (itoy%(int(0.1*toys.size()))==0)
       cout<<" Chebychev fit toy "<<itoy<<" / "<<toys.size()<<endl;
     TString stoy(to_string(itoy));
     RooRealVar nBkg_cheb_toy("nbkg_toy_"+stoy+name,"",sumEntries,0,2*sumEntries);
     RooRealVar nSgn_cheb_toy("nsgn_toy_"+stoy+name,"",0,-1*range_nSgn,range_nSgn);
     RooChebychev * cheb_pdf_toy = CreateChebychev( "cheb_toy_"+stoy+name, order, dilep_mass_out);    
     RooAddPdf total_pdf_cheb_toy("total_pdf_cheb_"+stoy,"", RooArgList(* cheb_pdf_toy, * sgn_pdf),  RooArgList(nBkg_cheb_toy,nSgn_cheb_toy));
     RooFitResult* toy_cheb_result = total_pdf_cheb_toy.fitTo(*toys[itoy],RooFit::Extended(1),RooFit::Save(),RooFit::Range(70 , 110),RooFit::PrintLevel(-1));
     hpull_cheb->Fill((nSgn_cheb_toy.getVal()-nExpected)/nSgn_cheb_toy.getError());
     if ( PrintPulls) 
        cout<<" Cheb3: #toy "<<itoy<<" pull: "<<(nSgn_cheb_toy.getVal()-nExpected)/nSgn_cheb_toy.getError()<<" Fit: nSgn "<<nSgn_cheb_toy.getVal()<<" +/- "<<nSgn_cheb_toy.getError()<<" nBkg "<<nBkg_cheb_toy.getVal()<<" +/- "<<nBkg_cheb_toy.getError()<<" entries "<<"; Gen Sgn "<<nExpected<<endl;
  }
  return hpull_cheb;
} 

#endif
