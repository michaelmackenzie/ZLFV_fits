#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// Z prime v2:
//////////////////// +based on v15
////////////////////



int ScanMuE_fit_sgn_v2(TString name="bin1_r2", float min_mass=100,
    float max_mass=500, float sgn_central=110, float sgn_width=2.08,
    float expected_Nsgn=65, bool create_dc_input=true, TString outvar="mass_ll",
    bool syst_sgn=false, TString varname="bin"){

   //////////////////////////////////// configuration /////////////////////////
   gROOT->SetBatch(true);
   TString outdir = "WorkspaceScanSGN/";
   gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", outdir.Data(), outdir.Data()));
   figdir_ = "figures/" + name + "/";
   gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir_.Data(), figdir_.Data()));

   bool GaussFnc=true;

   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ScanMuE fit v2: signal part ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;

   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", (max_mass+min_mass)/2., min_mass , max_mass, "GeV/c^{2}");

   RooRealVar mean("mean_"+varname,"",sgn_central);
   RooRealVar width("width_"+varname,"",sgn_width);
   mean.setConstant(true);
   width.setConstant(true);

   const bool use_energy_scale = true; //include energy scale uncertainties
   RooRealVar elec_ES_shift("elec_ES_shift_"+varname, "Electron ES shift", 0., -5., 5.); elec_ES_shift.setConstant(true);
   RooRealVar muon_ES_shift("muon_ES_shift_"+varname, "Muon ES shift"    , 0., -5., 5.); muon_ES_shift.setConstant(true);
   //Define the electron scale uncertainty as 0.5% and muon as 0.3%, assume the effect ~sqrt(sys/nom)
   // --> mean shift = 0.25% for electron ES and 0.15% for muon ES
   RooFormulaVar mean_func("mean_func_"+varname, "mean with offset",
                           "@0*(1 + @1*0.00250 + @2*0.00150)", RooArgList(mean, elec_ES_shift, muon_ES_shift));

   RooAbsPdf * signal_pdf;

   if (GaussFnc)
     signal_pdf = new RooGaussian("signal_pdf_"+varname,"gauss",dilep_mass_out,mean_func,width);

   RooRealVar n_sgn("signal_pdf_"+varname+"_norm","",expected_Nsgn,0,100000);
   n_sgn.setConstant(true);

   RooAddPdf signal_epdf("signal_epdf_"+varname,"ext_gauss",RooArgList(*signal_pdf),RooArgList(n_sgn));

   auto sgn_frame = dilep_mass_out.frame();
   signal_epdf.plotOn(sgn_frame,RooFit::LineColor(kBlue),RooFit::Name("sgn_epdf"),RooFit::Normalization(n_sgn.getVal(), RooAbsReal::NumEvent));
   save_plot(sgn_frame,"m(#mu,e)","scansgn_v2_gauss_"+name);

   //free the energy scale uncertainties if being included
   if(use_energy_scale) {
     elec_ES_shift.setConstant(false);
     muon_ES_shift.setConstant(false);
   }

   if (create_dc_input){
       RooWorkspace *wspace_sgn = new RooWorkspace("workspace_signal","workspace_signal");
       wspace_sgn->import(*signal_pdf);
       wspace_sgn->import(n_sgn);
       wspace_sgn->writeToFile(outdir+"workspace_scansgn_v2_"+name+".root");
    }

 return 0;
}
