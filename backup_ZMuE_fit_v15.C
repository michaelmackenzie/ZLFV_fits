#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// v15 updates:
//////////////////// +add more options for data obs 
//////////////////// +clean code and introduce fit_helper_core
//////////////////// +fix ndof
//////////////////// +ftest implemented 
//////////////////// +external functions for toys
//////////////////// +add th1 as data option



int ZMuE_fit_v15(TString name="bin1_r2", 
  //  TString data_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root",
    TString data_file="pseudo_data_from_MC_v2_r0.root",
    TString sgn_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root",
    TString xgbmin="0.3",TString xgbmax="0.7", bool plot_primitives=true,
    bool unblind=false, bool create_dc_input=false, float expected_Nsgn=0,
    TString outvar="mass_ll", bool syst_sgn=false, bool altfit_bkg=true, 
    TString varname="bin", bool pseudodata_input=false, 
    float pseudodata_norm=-1.0, int ntoys=-1, bool histo_input=false){

   //////////////////////////////////// configuration /////////////////////////
   gROOT->SetBatch(true);
   TString sgn_mu_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEUP_bdt_v7_signal_mcRun1*.root";
   TString sgn_mu_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEDOWN_bdt_v7_signal_mcRun1*.root"; 
   TString sgn_ele_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEUP_bdt_v7_signal_mcRun1*.root";
   TString sgn_ele_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEDOWN_bdt_v7_signal_mcRun1*.root";
   TString dilep_var_name="mass_ll";
   TString dilep_muup_name="mass_ll_Muon_scale_up";
   TString dilep_mudown_name="mass_ll_Muon_scale_down";
   TString dilep_eleup_name="mass_ll_Electron_scale_up";
   TString dilep_eledown_name="mass_ll_Electron_scale_down";
   std::vector<TString> syst_names{"Muon_scale_up", "Muon_scale_down", "Electron_scale_up", "Electron_scale_down"};
   std::vector<TString> syst_files{sgn_mu_up_file, sgn_mu_down_file, sgn_ele_up_file, sgn_ele_down_file};
  

  TString data_combine_file="template_zemu_v2_custom_ZmmF_1.05"; //IF histogram is used as  template (option histo_template_combine):  last part of the file MUST be the _"varname" option. Do not add the full file name juat up to _"varname"

   bool pseudodata_template_combine=false;
   bool pseudodata_fit_combine=false;
   bool data_blinded_combine=false;
   bool histo_template_combine=true; 
        TString histo_name_for_template="hbkg_safe"; // histo name

   TString bkg_pull_file="template_zemu_embed_v2"; //can be either MC ttree of blg or histo of bkg. IF histogram is used as  template (option histo_template_combine):  last part of the file MUST be the _"varname" option. Do not add the full file name juat up to _"varname"
   bool histo_template_pull=true;
        TString histo_name_for_template_pull="hbkg_safe"; // histo name
   float signal_toy_r=0.0;  
   bool Save_fittoy_plots=false;
   bool Print_fittoy_pull=false;
 
   int min_cheb_order=3,max_cheb_order=5;
   bool altfit_bkg_bst3=false;
   bool altfit_bkg_bst4=false;
   bool altfit_bkg_exp=false;
      int min_sumexp_order=1,max_sumexp_order=5;
      bool sumexp_recurse_coef=true;
   bool altfit_bkg_plaw=false;
      int min_sumplaw_order=1,max_sumplaw_order=5;
      bool sumplaw_recurse_coef=true;
 
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

   int nbin_blind =0;
   if ( !unblind )
      nbin_blind = (blind_max-blind_min) /( max_fit_range - min_fit_range )*nbin_data;

   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ZMuE_fit: v15 ***Main bkg: Chebychev3 ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;
   
   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   
   // read trees
   TTree * sgn_tree = get_tree("mytreefit",sgn_file,cuts);
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooDataSet sgn_dataset("sgn_dataset","sgn_dataset",RooArgSet(dilep_mass),RooFit::Import(*sgn_tree));
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");

   dilep_mass.setRange("left",min_fit_range, blind_min); //70,88
   dilep_mass.setRange("right",blind_max, max_fit_range); //94,110
   dilep_mass.setRange("sr",blind_min,blind_max); //88,94
   dilep_mass.setRange("full",min_fit_range, max_fit_range);


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
   RooRealVar n_sgn("sgn_PDF_norm", "",sgn_dataset.sumEntries(),-10*sgn_dataset.sumEntries(),10*sgn_dataset.sumEntries());
   RooAddPdf esgn_PDF("sgn","esgnPDF",RooArgList(sgn_PDF),RooArgList(n_sgn));

   ///// fit
   RooFitResult * sgn_result = esgn_PDF.fitTo(sgn_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range, max_fit_range),RooFit::PrintLevel(-1));
   ///// fix shape
   sgn_cb_mean.setConstant(true);   sgn_cb_width.setConstant(true);
   sgn_cb_a1.setConstant(true);     sgn_cb_n1.setConstant(true);
   sgn_cb_a2.setConstant(true);     sgn_cb_n2.setConstant(true);
   RooDoubleCB sgn_PDF_out("sgn_PDF_"+varname,"cb",dilep_mass_out,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);

   auto sgn_frame = dilep_mass.frame();   
   print_details(sgn_result);
   int n_param_sgn = sgn_result->floatParsFinal().getSize();
   float chi2_sgn = get_chi_squared(dilep_mass, &esgn_PDF, sgn_dataset, true, nbin_data, n_param_sgn);
   sgn_dataset.plotOn(sgn_frame,RooFit::Binning(nbin_data));
   esgn_PDF.plotOn(sgn_frame,RooFit::LineColor(kBlue));

   n_sgn.setVal(expected_Nsgn);
   if (plot_primitives)
      save_plot(sgn_frame,"m(#mu,e)","v15_prmtv_sgn_"+name);
   
   cout<<" Result  whole range    |    85-95 only"<<endl;
   std::pair<double,double> nSgn = yield_calc( n_sgn.getVal(), dilep_mass, &sgn_PDF);
   cout<<" - Expected nSgn "<<nSgn.first<<"   |  "<<nSgn.second<<endl;
   cout<<" ************************************************* "<<endl;
    
   /////// signal systematics
   double mean_shift_mu_up=0,width_shift_mu_up=0,mean_shift_mu_down=0,width_shift_mu_down=0;
   double mean_shift_ele_up=0,width_shift_ele_up=0,mean_shift_ele_down=0,width_shift_ele_down=0;
   
   if (syst_sgn){
     for (int isyst=0; isyst<syst_files.size(); isyst++){
        std::pair<double,double> max_mean_width = SignalSystematicsMaxMeanWidth(syst_files[isyst], cuts, syst_names[isyst], min_fit_range, max_fit_range, nbin_data, name) ;
        cout<<"signal systematic "+syst_names[isyst]<<" mean "<<max_mean_width.first<<" width "<<max_mean_width.second<<endl;
        if (isyst==0){ mean_shift_mu_up=max_mean_width.first; width_shift_mu_up=max_mean_width.second;}
        if (isyst==1){ mean_shift_mu_down=max_mean_width.first; width_shift_mu_down=max_mean_width.second;}
        if (isyst==2){ mean_shift_ele_up=max_mean_width.first; width_shift_ele_up=max_mean_width.second;}
        if (isyst==3){ mean_shift_ele_down=max_mean_width.first; width_shift_ele_down=max_mean_width.second;}
     }
   }
   /////// create signal workspace
   if (create_dc_input){
      RooRealVar n_sgn_out("sgn_PDF_"+varname+"_norm", "",expected_Nsgn,-10*expected_Nsgn,10*expected_Nsgn);
      n_sgn_out.setConstant(true);
     
      RooWorkspace *wspace_sgn = new RooWorkspace("workspace_signal","workspace_signal");
      wspace_sgn->import(sgn_PDF_out);
      wspace_sgn->import(n_sgn_out);

      double max_shift_muon=0;
      double max_shift_electron=0;
      if (syst_sgn){
        max_shift_muon= (fabs(sgn_cb_mean.getVal()-mean_shift_mu_down) + fabs(sgn_cb_mean.getVal()-mean_shift_mu_up) )/(2*sgn_cb_mean.getVal());
        max_shift_electron= fabs(mean_shift_ele_down-mean_shift_ele_up)/(2*sgn_cb_mean.getVal());
        cout<<"delta/M muon "<<max_shift_muon<<" electron "<<max_shift_electron<<endl;  
      }
      RooRealVar *mu_scale = new RooRealVar("nuisance_mu_scale_"+varname,"",0,-1,1);
      RooRealVar *ele_scale = new RooRealVar("nuisance_ele_scale_"+varname,"",0,-1,1);
      mu_scale->setConstant(true);
      ele_scale->setConstant(true);
      RooFormulaVar *mean_formula=new RooFormulaVar("mean_formula_"+varname, "", "@0*(1+"+TString(std::to_string(max_shift_muon))+"*@1)*(1+"+TString(std::to_string(max_shift_muon))+"*3*@2)", RooArgList(sgn_cb_mean,*mu_scale,*ele_scale));
      RooDoubleCB syst_sgn_PDF_out ("syst_sgn_PDF_"+varname,"cb",dilep_mass_out,*mean_formula,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
      RooRealVar n_syst_sgn_out("syst_sgn_PDF_"+varname+"_norm", "",expected_Nsgn,-10*expected_Nsgn,10*expected_Nsgn);
      n_syst_sgn_out.setConstant(true);
      if (syst_sgn){
        wspace_sgn->import(syst_sgn_PDF_out);
        wspace_sgn->import(n_syst_sgn_out);
      }
      wspace_sgn->writeToFile("workspace_v15_sgn_"+name+".root");
   }
   
   ///////////////////////////////// Data/bkg fit ////////////////////////////
   cout<<"\n *********************** Data fit ********************** "<<endl;

   RooRealVar NormGen_wt("NormGen_wt","gen weight", 0,-100,100);
   RooFormulaVar wgtFunc("wgtFunc","weight formula","NormGen_wt*"+TString(std::to_string(pseudodata_norm)),NormGen_wt);
   RooDataSet dataset = GetDataSet(data_file, cuts, dilep_mass, pseudodata_input, NormGen_wt,wgtFunc,pseudodata_norm,histo_input);

   cout<<"dataset = "<<dataset.sumEntries()<<endl;
   if (dataset.isWeighted()) cout<<"Is weighted"<<endl;
   else cout<<"Is NOT weighted"<<endl;
   
   ws.import(dataset); 

   std::vector<RooAbsPdf*> bkg_pdfs;
   std::vector<RooRealVar*> bkg_ampl;
   std::vector<TString> bkg_fnc_names;
   std::vector<TString> bkg_fnc_legs;


   /////// chebychev (main)
   //// define
   std::vector<RooAbsPdf*> bkg_cheb_pdfs;
   std::vector<RooRealVar*> bkg_cheb_ampl;
   std::vector<TString> bkg_cheb_names;
   std::vector<TString> bkg_cheb_legs;
   std::vector<int> bkg_cheb_orders;
   for (int iorder=min_cheb_order; iorder<max_cheb_order+1; iorder++){
      TString sorder(std::to_string(iorder));
      bkg_cheb_pdfs.push_back(CreateChebychev( "cheb"+sorder+"_bkgPDF", iorder, dilep_mass));
      bkg_cheb_ampl.push_back(new RooRealVar("cheb"+sorder+"_bkgPDF_norm","",dataset.sumEntries(),0,2*dataset.sumEntries()));
      bkg_cheb_names.push_back("cheb"+sorder+"_"+name);
      bkg_cheb_legs.push_back("Chebychev "+sorder);
      bkg_cheb_orders.push_back(iorder);
   }
   ///// ftest
   cout<<" ************************ Ftest Cheby begin ************************ "<<endl;
   std::pair<int,float> cheb_best_order =  Ftest(bkg_cheb_pdfs, bkg_cheb_ampl, dataset,dilep_mass, bkg_cheb_orders, bkg_cheb_names, bkg_cheb_legs, nbin_data, nbin_blind, (pseudodata_input || (!pseudodata_input && unblind)), pseudodata_input, name+"_cheb", 0.05,0.01,printout_levels);
   cout<<" ************************ Ftest Cheby end ************************ "<<endl;
   ///// best chebychev
   bkg_pdfs.push_back(CreateChebychev( "cheb_bkgPDF", cheb_best_order.first, dilep_mass));
   bkg_ampl.push_back(new RooRealVar("cheb_bkgPDF_norm","",dataset.sumEntries(),0,2*dataset.sumEntries()));
   bkg_fnc_names.push_back("cheb_"+name);
   bkg_fnc_legs.push_back("Chebychev "+TString(std::to_string( cheb_best_order.first)));
   


   /////// 3rd order B. 
   std::vector<RooRealVar> bkg_bst3_params= BstParams(3, varname);
   RooBernsteinFast<3> bst3_bkgPDF("bst3_bkgPDF_"+varname, "", dilep_mass, RooArgList(bkg_bst3_params[0],bkg_bst3_params[1],bkg_bst3_params[2]) );
   bst3_bkgPDF.protectSubRange(true);
   RooRealVar bst3_n_bkg("bst3_bkgPDF_norm","",dataset.sumEntries(),0,2*dataset.sumEntries());   
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
   RooRealVar bst4_n_bkg("bst4_bkgPDF_norm","",dataset.sumEntries(),0,2*dataset.sumEntries());   
   if (altfit_bkg && altfit_bkg_bst4){
      bkg_pdfs.push_back(&bst4_bkgPDF);
      bkg_ampl.push_back(&bst4_n_bkg);
      bkg_fnc_names.push_back("bst4_"+name);
      bkg_fnc_legs.push_back("Bernstein 4");
    }

   /////// sum exp 
   std::pair<int,float> sumexp_best_order(0,0);
   if (altfit_bkg && altfit_bkg_exp){
     ///// define
     std::vector<RooAbsPdf*> bkg_sumexp_pdfs;
     std::vector<RooRealVar*> bkg_sumexp_ampl;
     std::vector<TString> bkg_sumexp_names;
     std::vector<TString> bkg_sumexp_legs;
     std::vector<int> bkg_sumexp_orders;
     for (int iorder=min_sumexp_order; iorder<max_sumexp_order+1; iorder++){
        TString sorder(std::to_string(iorder));
        bkg_sumexp_pdfs.push_back(CreateSumExpo( "sumexp"+sorder+"_bkgPDF", iorder, dilep_mass, sumexp_recurse_coef));
        bkg_sumexp_ampl.push_back(new RooRealVar("sumexp"+sorder+"_bkgPDF_norm","",dataset.sumEntries(),0.5*dataset.sumEntries(),1.5*dataset.sumEntries()));
        bkg_sumexp_names.push_back("sumexp"+sorder+"_"+name);
        bkg_sumexp_legs.push_back("#Sigma Expo "+sorder);
        bkg_sumexp_orders.push_back(iorder);
     }
     ////// ftest
     cout<<" ************************ Ftest SumExp begin ************************ "<<endl;
     sumexp_best_order =  Ftest(bkg_sumexp_pdfs, bkg_sumexp_ampl, dataset,dilep_mass, bkg_sumexp_orders, bkg_sumexp_names, bkg_sumexp_legs, nbin_data, nbin_blind-1, (pseudodata_input || (!pseudodata_input && unblind)), pseudodata_input, name+"_sumexp", 0.05,0.01,printout_levels);
     cout<<" ************************ Ftest SumExp end ************************ "<<endl;
     ////// best sumexp
     bkg_pdfs.push_back(CreateSumExpo( "sumexp_bkgPDF", sumexp_best_order.first, dilep_mass,sumexp_recurse_coef));
     bkg_ampl.push_back(new RooRealVar("sumexp_bkgPDF_norm","",dataset.sumEntries(),0,2*dataset.sumEntries()));
     bkg_fnc_names.push_back("sumexpo_"+name);
     bkg_fnc_legs.push_back("#Sigma expo "+TString(to_string(sumexp_best_order.first)));
   }


   /////// sum power law
   std::pair<int,float> sumplaw_best_order(0,0);
   if (altfit_bkg && altfit_bkg_plaw){
     ///// define
     std::vector<RooAbsPdf*> bkg_sumplaw_pdfs;
     std::vector<RooRealVar*> bkg_sumplaw_ampl;
     std::vector<TString> bkg_sumplaw_names;
     std::vector<TString> bkg_sumplaw_legs;
     std::vector<int> bkg_sumplaw_orders;
     for (int iorder=min_sumplaw_order; iorder<max_sumplaw_order+1; iorder++){
        TString sorder(std::to_string(iorder));
        bkg_sumplaw_pdfs.push_back(CreateSumPower( "sumplaw"+sorder+"_bkgPDF", iorder, dilep_mass));
        bkg_sumplaw_ampl.push_back(new RooRealVar("sumplaw"+sorder+"_bkgPDF_norm","",dataset.sumEntries(),0.5*dataset.sumEntries(),2*dataset.sumEntries()));
        bkg_sumplaw_names.push_back("sumplaw"+sorder+"_"+name);
        bkg_sumplaw_legs.push_back("#Sigma power law "+sorder);
        bkg_sumplaw_orders.push_back(iorder);
     }
     ////// ftest
     cout<<" ************************ Ftest Sum Power Law begin ************************ "<<endl;
     sumplaw_best_order =  Ftest(bkg_sumplaw_pdfs, bkg_sumplaw_ampl, dataset, dilep_mass, bkg_sumplaw_orders, bkg_sumplaw_names, bkg_sumplaw_legs, nbin_data, nbin_blind, (pseudodata_input || (!pseudodata_input && unblind)), pseudodata_input, name+"_sumplaw", 0.05,0.01,printout_levels);
     cout<<" ************************ Ftest Sum Power Law end ************************ "<<endl;
     ////// best sumplaw
     bkg_pdfs.push_back(CreateSumPower( "sumplaw_bkgPDF", sumplaw_best_order.first, dilep_mass));
     bkg_ampl.push_back(new RooRealVar("sumplaw_bkgPDF_norm","",dataset.sumEntries(),0,2*dataset.sumEntries()));
     bkg_fnc_names.push_back("sumplaw_"+name);
     bkg_fnc_legs.push_back("#Sigma power law "+TString(to_string(sumplaw_best_order.first)));
   }



   /////////////////////////////// Fits ///////////////////////////////////
   //////// Bkg only
   std::vector<std::vector<float>> chi2_bkgfit = FitBkgFunctions(bkg_pdfs,bkg_ampl,dataset,dilep_mass,bkg_fnc_names,bkg_fnc_legs,Bkg_only_fit_whole_region, unblind, nbin_data,nbin_blind, pseudodata_input, name);

   std::vector<float> cheb_best_parameters;
   for (int i=2; i<chi2_bkgfit[0].size(); i++)
     cheb_best_parameters.push_back(chi2_bkgfit[0][i]);
    
   std::vector<float> sumexp_best_parameters;
   std::vector<float> sumplaw_best_parameters;
   if (altfit_bkg){ 
     int nfamily=1;
     if (altfit_bkg_exp) {
       for (int i=2; i<chi2_bkgfit[nfamily].size(); i++)
         sumexp_best_parameters.push_back(chi2_bkgfit[nfamily][i]);
       nfamily+=1;
     }    
     if (altfit_bkg_plaw){
       for (int i=2; i<chi2_bkgfit[nfamily].size(); i++)
         sumplaw_best_parameters.push_back(chi2_bkgfit[nfamily][i]);
     }
   }


   ///////// unblinding
   if (unblind)
     std::vector<std::pair<float,int>> chi2_totalfit = FitTotalFunctions(bkg_pdfs, &sgn_PDF, dilep_mass, dataset, bkg_fnc_names, nbin_data, pseudodata_input, name);
   

    
   //////////////////// create bkg and data workspace ////////////////////////
   if (create_dc_input){  
     dilep_mass_out.setBins(nbin_data);
     cout<<" ******************* creating data card ****************** "<<endl;
     RooDataSet * data_obs = GetDataObs(dilep_mass_out, pseudodata_fit_combine, data_blinded_combine, histo_template_combine, pseudodata_template_combine, data_combine_file, cuts, data_combine_file+"_"+varname+".root", histo_name_for_template, min_fit_range, max_fit_range, blind_min, blind_max, dataset.sumEntries(), "v15", name);

     data_obs->SetName("sim_data_obs_"+varname);
     RooDataHist hdata_sim("sim_binned_obs_"+varname,"sim_binned_obs_"+varname,dilep_mass_out,*data_obs); // create also binned version

     auto dobs_frame = dilep_mass_out.frame();
     data_obs->plotOn(dobs_frame,RooFit::MarkerColor(1),RooFit::LineColor(1));
     save_plot(dobs_frame,"m(e,#mu)","v15_data_obs_"+name);

     /// generate signal toys
     RooDataSet * signal_obs_r1 = sgn_PDF_out.generate(RooArgSet(dilep_mass_out),expected_Nsgn);
     RooDataSet *data_obs_r1 =(RooDataSet*) data_obs->Clone("sim_data_obs_r1_"+varname);
     data_obs_r1->append( *signal_obs_r1);

     ///// output bkg functions
     RooCategory cat("pdfindex_"+varname, "");
     RooArgList models_out; //container for multpdf
     RooWorkspace *wspace_bkg = new RooWorkspace("workspace_background","workspace_background"); //background & data workspace
     RooWorkspace *wspace_singles_bkg = new RooWorkspace("workspace_singles_background","workspace_singles_background"); //background / data workspace

     // cheb output
     RooChebychev* cheb_bkgPDF_out = CreateChebychev( "cheb"+TString(to_string(cheb_best_order.first))+"_bkgPDF_"+varname,cheb_best_order.first , dilep_mass_out,cheb_best_parameters);
     RooRealVar cheb_n_bkg_out("cheb"+TString(to_string(cheb_best_order.first))+"_bkgPDF_"+varname+"_norm","",dataset.sumEntries(),0,2*dataset.sumEntries());
     models_out.add(*cheb_bkgPDF_out);
     wspace_singles_bkg->import(*cheb_bkgPDF_out);
     wspace_singles_bkg->import(cheb_n_bkg_out);

     // sum exp output
     RooAddPdf * sumexp_bkgPDF_out = CreateSumExpo("sumexp"+TString(to_string(sumexp_best_order.first))+"_bkgPDF_"+varname, sumexp_best_order.first, dilep_mass_out,sumexp_recurse_coef,sumexp_best_parameters);
     RooRealVar sumexp_n_bkg_out("sumexp"+TString(to_string(sumexp_best_order.first))+"_bkgPDF_"+varname+"_norm","",dataset.sumEntries(),0,2*dataset.sumEntries());
     if (altfit_bkg && altfit_bkg_exp){
        models_out.add(*sumexp_bkgPDF_out);
        wspace_singles_bkg->import(*sumexp_bkgPDF_out);
	wspace_singles_bkg->import(sumexp_n_bkg_out);
     }

     // sum plaw output
     RooAddPdf * sumplaw_bkgPDF_out = CreateSumPower("sumplaw"+TString(to_string(sumplaw_best_order.first))+"_bkgPDF_"+varname, sumplaw_best_order.first, dilep_mass_out,sumplaw_best_parameters);
     RooRealVar sumplaw_n_bkg_out("sumplaw"+TString(to_string(sumplaw_best_order.first))+"_bkgPDF_"+varname+"_norm","",dataset.sumEntries(),0,2*dataset.sumEntries());
     if (altfit_bkg && altfit_bkg_plaw){
        models_out.add(*sumplaw_bkgPDF_out);
        wspace_singles_bkg->import(*sumplaw_bkgPDF_out);
	wspace_singles_bkg->import(sumplaw_n_bkg_out);
     }


     // import pdfs to wspace
     RooMultiPdf multipdf("multipdf_"+varname, "", cat, models_out);
//     multipdf.setCorrectionFactor(0.0001);
     RooRealVar norm_out("multipdf_"+varname+"_norm","",dataset.sumEntries(),0.5*dataset.sumEntries(),1.5*dataset.sumEntries());  
    
     wspace_bkg->import(cat);
     wspace_bkg->import(multipdf);
     wspace_bkg->import(norm_out);
   
     wspace_bkg->import(*data_obs);
     wspace_bkg->import(*data_obs_r1);
     wspace_bkg->import(hdata_sim);

     // mc template as additional pdf
     if ( pseudodata_fit_combine || pseudodata_template_combine || histo_template_combine ){
        TH1F* htemplate = GetHistoTemplate(data_combine_file, cuts, min_fit_range, max_fit_range, histo_name_for_template, data_combine_file+"_"+varname+".root", histo_template_combine,"pdf" );

        RooDataHist rdh_template("rdh_template_"+varname,"",RooArgSet(dilep_mass_out),htemplate);
        RooHistPdf pdf_template("pdf_template_"+varname,"",RooArgSet(dilep_mass_out),rdh_template);
        RooRealVar pdf_template_norm("pdf_template_"+varname+"_norm","",dataset.sumEntries(),0,2*dataset.sumEntries());
        wspace_singles_bkg->import(pdf_template);
        wspace_singles_bkg->import(pdf_template_norm);
     }
     wspace_bkg->writeToFile("workspace_v15_bkg_"+name+".root"); // write outputt
     wspace_singles_bkg->writeToFile("workspace_v15_bkg_single_fnc_"+name+".root"); 
   }
    
    
   ///////////////////////////////////  toy tests /////////////////////////////
   if (ntoys>0){

     cout<<" ****************  starting pull tests ******************"<<endl;

     TH1F* htemplate_pull = GetHistoTemplate(bkg_pull_file, cuts, min_fit_range, max_fit_range, histo_name_for_template_pull, bkg_pull_file+"_"+varname+".root", histo_template_pull,"pull" );
     
     save_th1({htemplate_pull}, {""},"m(e#mu)","v15_pull_template");
     RooDataHist rd_htoy("rd_htoy","rd_htoy",RooArgSet(dilep_mass_out),htemplate_pull);
     RooHistPdf *pdf_htoy = new RooHistPdf("pdf_htoy","pdf_htoy",RooArgSet(dilep_mass_out),rd_htoy);
     
     std::vector<RooDataSet *> toy_datasets = GenerateTemplateToys( pdf_htoy, (&sgn_PDF_out), dilep_mass_out, ntoys, dataset.sumEntries(), expected_Nsgn*signal_toy_r);

     TH1F* hpull_cheb = ChebychevPullFromToys( cheb_best_order.first, (&sgn_PDF_out), dilep_mass_out, toy_datasets, expected_Nsgn*signal_toy_r, dataset.sumEntries(), "order"+TString(to_string(cheb_best_order.first)),Print_fittoy_pull);

     save_th1({hpull_cheb}, {""}, "pull", "v15_pull_"+name+"_toyR"+TString(to_string(signal_toy_r)), true);

   }

 return 0;
} 
