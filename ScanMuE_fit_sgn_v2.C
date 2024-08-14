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

   bool GaussFnc=true;

   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ScanMuE fit v2: signal part ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;

   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_mass , max_mass, "GeV/c^{2}");

   RooRealVar mean("mean_"+varname,"",sgn_central);
   RooRealVar width("width_"+varname,"",sgn_width);
   mean.setConstant(true);
   width.setConstant(true);
   
   RooAbsPdf * signal_pdf;

   if (GaussFnc)
     signal_pdf = new RooGaussian("signal_pdf_"+varname,"gauss",dilep_mass_out,mean,width);

   RooRealVar n_sgn("signal_pdf_"+varname+"_norm","",expected_Nsgn,0,100000);
   n_sgn.setConstant(true);

   RooAddPdf signal_epdf("signal_epdf_"+varname,"ext_gauss",RooArgList(*signal_pdf),RooArgList(n_sgn));

   auto sgn_frame = dilep_mass_out.frame();
   signal_epdf.plotOn(sgn_frame,RooFit::LineColor(kBlue),RooFit::Name("sgn_epdf"),RooFit::Normalization(n_sgn.getVal(), RooAbsReal::NumEvent));
   save_plot(sgn_frame,"m(#mu,e)","WorkspaceScanSGN/scansgn_v2_gauss_"+name);
   
   if (create_dc_input){
       RooWorkspace *wspace_sgn = new RooWorkspace("workspace_signal","workspace_signal");
       wspace_sgn->import(*signal_pdf);
       wspace_sgn->import(n_sgn);
       wspace_sgn->writeToFile("WorkspaceScanSGN//workspace_scansgn_v2_"+name+".root");
    }   
		    
 return 0;
} 
