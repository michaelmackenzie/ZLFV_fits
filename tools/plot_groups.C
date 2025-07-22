//Plot uncertainty impacts from standard groups
#include "extract_impact.C"
int verbose_ = 1;

int plot_groups(const char* file, const char* tag = "_zmutau_v09j", const bool doObs = false, const bool run = true,
                TString only_group = "") {

  vector<TString> groups = {
    "elec_ES_shift",
    "muon_ES_shift",
    "BTag",
    "Pileup",
    "Lumi",
    "Theory",
    "g_JetMET",
    "g_IDs",
    "g_MCStats",
    "g_Env",
    "g_PDF",
    "All_Systematics"
  };

  if(run) {
    TString args = "-M FitDiagnostics --rMin -50 --rMax 50 --cminDefaultMinimizerStrategy 0";
    args += " --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants";
    args += " --cminRunAllDiscreteCombinations --X-rtd MINIMIZER_multiMin_maskChannels=2";
    args += " --cminApproxPreFitTolerance 0.01 --cminPreScan --cminPreFit 1";
    args += " --cminDefaultMinimizerTolerance 0.001 --cminDiscreteMinTol 0.0001";
    cout << "\n>>> Performing the nominal fit\n";
    TString command = Form("combine %s -n _groupFit_Test%s_Nominal  -d %s %s", args.Data(), tag, file, (doObs) ? "" : "-t -1");
    if(only_group == "") {
      if(verbose_) cout << command.Data() << endl;
      gSystem->Exec(command.Data());
    }
    for(TString group : groups) {
      if(only_group != "" && group != only_group) continue;
      const char* output_name = Form("fitDiagnostics_groupFit_Test%s_%s.root", tag, group.Data());
      //remove previous fit results if there
      gSystem->Exec(Form("[ -f %s ] && rm %s", output_name, output_name));
      if(group == "All_Systematics") {
        cout << "\n>>> Fitting without any systematics\n";
        TString command = Form("combine %s -n _groupFit_Test%s_%s  -d %s %s --freezeParameters allConstrainedNuisances",
                               args.Data(), tag, group.Data(), file, (doObs) ? "" : "-t -1");
        if(verbose_) cout << command.Data() << endl;
        gSystem->Exec(command.Data());
      } else if(group == "g_Env") { // handle this individually
        cout << "\n>>> Fitting without group " << group.Data() << endl;
        TString command = Form("combine %s -n _groupFit_Test%s_%s  -d %s %s --freezeParameters pdfindex_bin1,pdfindex_bin2",
                               args.Data(), tag, group.Data(), file, (doObs) ? "" : "-t -1");
        if(verbose_) cout << command.Data() << endl;
        gSystem->Exec(command.Data());
      } else if(group == "g_PDF") { // handle this individually
        cout << "\n>>> Fitting without group " << group.Data() << endl;
        TString command = Form("combine %s -n _groupFit_Test%s_%s  -d %s %s --freezeParameters var{bkg_.*},pdfindex_bin1,pdfindex_bin2",
                               args.Data(), tag, group.Data(), file, (doObs) ? "" : "-t -1");
        if(verbose_) cout << command.Data() << endl;
        gSystem->Exec(command.Data());
      } else if(group.BeginsWith("g_")) {
        cout << "\n>>> Fitting without group " << group.Data() << endl;
        TString command = Form("combine %s -n _groupFit_Test%s_%s  -d %s %s --freezeNuisanceGroups %s",
                               args.Data(), tag, group.Data(), file, (doObs) ? "" : "-t -1", group.Data());
        if(verbose_) cout << command.Data() << endl;
        gSystem->Exec(command.Data());
      } else {
        cout << "\n>>> Fitting without parameter " << group.Data() << endl;
        TString command = Form("combine %s -n _groupFit_Test%s_%s  -d %s %s --freezeParameters %s",
                               args.Data(), tag, group.Data(), file, (doObs) ? "" : "-t -1", group.Data());
        if(verbose_) cout << command.Data() << endl;
        gSystem->Exec(command.Data());
      }
    }
  }

  //Add a pure statistical uncertainty group, approximated using the nominal and no systematics evaluations
  groups.push_back("Statistical");


  const int n = groups.size();
  double rs[n+1], ups[n+1], downs[n+1];
  double ys[n+1], yup[n+1], ydown[n+1];
  ys   [n] = n+1; //total fit results
  ydown[n] = 0.1;
  yup  [n] = 0.1;
  fiteffect_t effect;
  double max_val(-1.e10), min_val(1.e10);
  for(int index = 0; index < n; ++index) {
    TString group = groups[index];
    cout << "Evaluating group " << group.Data() << endl;
    ys   [index] = index+1;
    ydown[index] = 0.1;
    yup  [index] = 0.1;

    if(group == "Statistical") {
      //Get information for nominal and all systematics impact
      const int nom_index = n;
      const int sys_index = n - 2; //FIXME: All_Systematics assumed second to last for now
      const double nom_up   = ups  [nom_index];
      const double nom_down = downs[nom_index];
      const double sys_up   = ups  [sys_index];
      const double sys_down = downs[sys_index];
      // const double nom_up   = ups  [nom_index] - rs[nom_index];
      // const double nom_down = downs[nom_index] - rs[nom_index];
      // const double sys_up   = ups  [sys_index] - rs[sys_index];
      // const double sys_down = downs[sys_index] - rs[sys_index];
      rs   [index] = rs[nom_index];
      ups  [index] = std::sqrt(max(0., std::pow(nom_up  ,2) - std::pow(sys_up  , 2)));
      downs[index] = std::sqrt(max(0., std::pow(nom_down,2) - std::pow(sys_down, 2)));
      printf("Estimating statistical as: %.4f -%.4f +%.4f\n", rs[index], ups[index], downs[index]);
    } else {
      int status = extract_impact(Form("fitDiagnostics_groupFit_Test%s_Nominal.root", tag), Form("fitDiagnostics_groupFit_Test%s_%s.root", tag, group.Data()), &effect);
      if(status) {
        rs   [index] = -999;
        ups  [index] = -999;
        downs[index] = -999;
        continue;
      }
      rs   [index] = effect.r_group_;
      ups  [index] = effect.impact_up_;
      downs[index] = effect.impact_down_;
      rs   [n]     = effect.r_nom_;
      ups  [n]     = effect.up_nom_;
      downs[n]     = effect.down_nom_;
    }
    max_val = max(max_val, rs[index]+ups[index]);
    max_val = max(max_val, rs[n]+ups[n]);
    min_val = min(min_val, rs[index]-downs[index]);
    min_val = min(min_val, rs[n]-downs[n]);
  }

  // Print summary info
  printf("Summary info        :   r     up    down  (up/tot down/tot)\n");

  for(int index = 0; index <= n; ++index) {
    TString group = (index < n) ? groups[index] : "Total";
    printf("%-20s: %.3f +%.3f -%.3f (%.3f -%.3f)\n", group.Data(), rs[index], ups[index], downs[index], ups[index]/ups[n], downs[index]/downs[n]);
  }


  TGraph* g = new TGraphAsymmErrors(n+1, rs, ys, ups, downs, yup, ydown);
  TCanvas c("c", "c", 700, 1000);
  c.SetTopMargin(0.1); c.SetBottomMargin(0.1);
  g->SetLineWidth(2);
  g->SetLineColor(kBlue);
  g->SetFillColor(kBlue);
  g->SetTitle("Systematic uncertainties; #sigma_{r};");
  g->Draw("APE2");
  const double xbuffer = 0.1*(max_val-min_val);
  g->GetXaxis()->SetRangeUser(min_val-xbuffer,max_val+xbuffer);
  g->GetYaxis()->SetRangeUser(0, n+2);
  g->GetYaxis()->SetLabelSize(0); //remove number labels


  //Draw the category labels
  TLatex label;
  label.SetNDC();
  label.SetTextFont(72);
  label.SetTextColor(kBlack);
  label.SetTextSize(0.05);
  label.SetTextAlign(13);
  label.SetTextAngle(0);
  label.SetTextSize(0.03);
  for(int index = 0; index <= n; ++index) {
    const double ystart = 0.1;
    const double yend   = 0.9;
    const double yloc = ystart + (yend-ystart)*(index+1.1)/(n+2);
    if(index < n) {
      TString name = groups[index];
      name.ReplaceAll("g_", "");
      label.DrawLatex(0.01, yloc, Form("%s", name.Data()));
    } else {
      label.DrawLatex(0.01, yloc, "Total");
    }
  }
  c.SetGridx();
  c.Modified(); c.Update();

  c.SaveAs(Form("groups%s.png", tag));
  return 0;
}
