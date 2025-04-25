//Plot B->e+mu fit results from fit diagnostics
//create a fit diagnostics root file via:
//$> combine -M FitDiagnostics -d <input card/workspace> --saveShapes --saveWithUncertainties [additional options]

bool unblind_      = false;
int  err_mode_     =  1   ; //errors in the pulls: 0: sqrt(data^2 + fit^2); 1: sqrt(data^2 - fit^2)
bool debug_        = false; //print debug info
bool do_single_    = false; //test printing a single histogram
bool do_sig_wt_    = false; //make a plot with S/(S+B) weighting
bool only_s_fit_   = true ; //only make the S+B fit plots
bool is_prelim_    = false; //is preliminary or not
TString file_type_ = "pdf";

//------------------------------------------------------------------------------------------
// Helper functions
double hmax(TH1* h) {
  double max_val = h->GetBinContent(1);
  for(int ibin = 2; ibin <= h->GetNbinsX(); ++ibin) max_val = max(max_val, h->GetBinContent(ibin));
  return max_val;
}

double hmin(TH1* h, double cutoff = 0.01) {
  double min_val = h->GetMaximum();
  for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    if(h->GetBinContent(ibin) < cutoff) continue;
    min_val = min(min_val, h->GetBinContent(ibin));
  }
  return max(cutoff, min_val);
}

double gmax(TGraph* g) {
  const int nbins = g->GetN();
  double max_val = 0.;
  for(int ibin = 0; ibin < nbins; ++ibin) {
    const double val = g->GetY()[ibin];
    max_val = max(max_val, val);
  }
  return max_val;
}

double gmin(TGraph* g, double cutoff = 0.01) {
  const int nbins = g->GetN();
  double min_val = -1;
  for(int ibin = 0; ibin < nbins; ++ibin) {
    const double val = g->GetY()[ibin];
    if(val < cutoff) continue;
    min_val = (min_val < 0.) ? val : min(min_val, val);
  }
  return max(cutoff, min_val);
}


//------------------------------------------------------------------------------------------
// Draw the CMS label
void draw_cms_label(float left_margin = 0.10) {
    //CMS prelim drawing
    TText cmslabel;
    cmslabel.SetNDC();
    cmslabel.SetTextColor(1);
    cmslabel.SetTextSize(0.09);
    cmslabel.SetTextAlign(11);
    cmslabel.SetTextAngle(0);
    cmslabel.SetTextFont(61);
    cmslabel.DrawText(left_margin + 0.04, 0.81, "CMS");
    if(is_prelim_) {
      cmslabel.SetTextFont(52);
      cmslabel.SetTextSize(0.76*cmslabel.GetTextSize());
      cmslabel.DrawText(left_margin + 0.04, 0.75, "Preliminary");
    }
}

void draw_luminosity(int year = -1) {
  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  // label.SetTextColor(1);
  label.SetTextSize(0.055);
  label.SetTextAlign(31);
  label.SetTextAngle(0);
  TString period = (year > 2000) ? Form("%i, ", year) : "";
  const double lum = (year == 2016) ? 36.33 : (year == 2017) ? 41.48 : (year == 2018) ? 59.83 : 137.64;
  label.DrawLatex(0.97, 0.915, Form("%s%.0f fb^{-1} (13 TeV)",period.Data(),lum));
}

void draw_category(TString tag, bool zprime, float left_margin = 0.10) {
  TString cat = "";
  if(zprime) {
    if(tag.Contains("bin_1")) cat = "0.3 < BDT < 0.7";
    else                      cat = "0.7 < BDT < 1.0";
  } else {
    if     (tag.Contains("bin_1")) cat = "0.3 < BDT < 0.7";
    else if(tag.Contains("bin_2")) cat = "0.7 < BDT < 0.9";
    else                           cat = "0.9 < BDT < 1.0";
  }
  TText label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.055);
  label.SetTextAlign(11);
  label.SetTextAngle(0);
  if(is_prelim_) label.DrawText(left_margin + 0.04, 0.69, cat.Data());
  else           label.DrawText(left_margin + 0.04, 0.75, cat.Data());
}

//------------------------------------------------------------------------------------------
// Scale the data
void scale_g(TGraphAsymmErrors* g, float scale = 1.f) {
  const int nbins = g->GetN();
  for(int ibin = 0; ibin < nbins; ++ibin) {
    g->SetPointY(ibin, g->GetPointY(ibin)*scale);
    g->SetPointEYhigh(ibin, g->GetEYhigh()[ibin]*scale);
    g->SetPointEYlow(ibin, g->GetEYlow()[ibin]*scale);
  }
}

//------------------------------------------------------------------------------------------
// Get a distribution from the directory list
TH1* get_hist(vector<TDirectoryFile*> dirs, const char* name, vector<double> weights = {}) {
  TH1* h = nullptr;
  int index = 0;
  for(auto dir : dirs) {
    TH1* h_tmp = (TH1*) dir->Get(name);
    if(!h_tmp) {
      // cout << __func__ << ": Histogram " << name << " not found in directory " << dir->GetName() << endl;
      return nullptr;
    }
    const double weight = (weights.size() > index) ? weights[index] : 1.f;
    if(!h) {
      h = (TH1*) h_tmp->Clone(Form("%s_Run2", name));
      h->Scale(weight);
    } else {
      h->Add(h_tmp, weight);
    }
    ++index;
  }
  return h;
}

//------------------------------------------------------------------------------------------
// Get the data distribution from the directory list
TGraphAsymmErrors* get_data(vector<TDirectoryFile*> dirs, vector<double> weights = {}) {
  TGraphAsymmErrors* g = nullptr;
  int index = 0;
  for(auto dir : dirs) {
    TGraphAsymmErrors* g_tmp = (TGraphAsymmErrors*) dir->Get("data");
    if(!g_tmp) {
      cout << __func__ << ": Data not found in directory " << dir->GetName() << endl;
      return nullptr;
    }
    const double weight = (weights.size() > index) ? weights[index] : 1.f;
    if(!g) {
      g = (TGraphAsymmErrors*) g_tmp->Clone("data_Run2");
      if(weights.size() > 0) { //scale the data
        for(int ipoint = 0; ipoint < g->GetN(); ++ipoint) {
          g->SetPointY(ipoint, weights[0]*g->GetPointY(ipoint));
          //use sqrt(N) for the errors
          g->SetPointEYhigh(ipoint, weight*sqrt(g->GetPointY(ipoint)));
          g->SetPointEYlow (ipoint, weight*sqrt(g->GetPointY(ipoint)));
        }

      }
    } else {
      //Add the data and errors
      for(int ipoint = 0; ipoint < g->GetN(); ++ipoint) {
        g->SetPointY(ipoint, g->GetPointY(ipoint) + weight*g_tmp->GetPointY(ipoint));
        //use sqrt(N) and add in quadrature
        g->SetPointEYhigh(ipoint, sqrt(pow(g->GetErrorYhigh(ipoint), 2) + weight*g_tmp->GetPointY(ipoint)));
        g->SetPointEYlow (ipoint, sqrt(pow(g->GetErrorYhigh(ipoint), 2) + weight*g_tmp->GetPointY(ipoint)));
      }
    }
    ++index;
  }
  if(g) {
    //no x errors
    for(int ipoint = 0; ipoint < g->GetN(); ++ipoint) {
      g->SetPointEXhigh(ipoint, 0.);
      g->SetPointEXlow (ipoint, 0.);
    }
  }

  return g;
}

//------------------------------------------------------------------------------------------
// Print an individual distribution
int print_hist(vector<TDirectoryFile*> dirs, TString tag, TString outdir, bool s_over_sb = false) {

  //get the weights for each region if requested
  vector<double> weights;
  if(s_over_sb && dirs.size() > 1) {
    double max_wt = 0.;
    for(auto dir : dirs) {
      auto hbkg = get_hist({dir}, "total_background");
      auto hsig = get_hist({dir}, "total_signal");
      const double peak_mass = hsig->GetBinCenter(hsig->GetMaximumBin());
      const double width = hsig->GetStdDev();
      const double nbkg = hbkg->Integral(hbkg->FindBin(peak_mass - 2.*width), hbkg->FindBin(peak_mass + 2.*width));
      const double nsig = fabs(hsig->Integral()); //use |S| in the ratio
      weights.push_back(max(1.e-5, nsig/(nbkg + nsig)));
      max_wt = max(max_wt, weights.back());
    }
    for(unsigned index = 0; index < weights.size(); ++index) {
      weights[index] *= 1./max_wt; //normalize to max of 1
      cout << "Category " << index << " weight = " << weights[index] << endl;
    }
  }

  //Get the fit results and the data
  TH1* hbackground         = get_hist(dirs, "total_background", weights);
  TH1* hbkg                = get_hist(dirs, "bkg", weights);
  TH1* hzmumu              = get_hist(dirs, "zmumu", weights);
  TH1* hzemu               = get_hist(dirs, "zemu", weights);
  TH1* hsignal             = get_hist(dirs, "total_signal", weights);
  TH1* htotal              = get_hist(dirs, "total", weights);
  TGraphAsymmErrors* gdata = get_data(dirs, weights);

  //Z->mumu names in ZLFV_fits datacards
  if(!hzmumu) hzmumu = get_hist(dirs, "zmm", weights);
  if(!hzemu ) hzemu  = get_hist(dirs, "sgn", weights);

  //Z prime scan names
  const bool zprime = !hbkg || !hzemu;
  if(!hbkg ) hbkg  = get_hist(dirs, "background", weights);
  if(!hzemu) hzemu = get_hist(dirs, "signal", weights);

  if(!hsignal || !hbackground || !htotal || !hzemu || !gdata) {
    cout << "Data not found for tag " << tag.Data() << endl;
    return 1;
  }

  //N(data) points, ensure it matches the background model
  const int nbins = gdata->GetN();
  if(nbins != hbackground->GetNbinsX()) {
    cout << "Data and background have different bin numbers!\n";
    return 2;
  }

  //Scale the inputs to remove the 1/bin width normalization
  const float scale = hbkg->GetBinWidth(1);
  scale_g(gdata,scale);
  hbackground->Scale(scale);
  hsignal->Scale(scale);
  hbkg->Scale(scale);
  hzemu->Scale(scale);
  htotal->Scale(scale);
  if(hzmumu) hzmumu->Scale(scale);

  //Determine the x-axis range to use
  const double xmin = htotal->GetBinLowEdge(htotal->FindFirstBinAbove(0.01));
  const double xmax = htotal->GetXaxis()->GetBinUpEdge(htotal->FindLastBinAbove(0.01));

  //Create the canvas to plot on
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas* c = new TCanvas("c", "c", 1000, 900);
  const float x1(0.23), x2(0.43);
  TPad* pad1 = new TPad("pad1", "pad1", 0., x2, 1., 1.);
  TPad* pad2 = new TPad("pad2", "pad2", 0., x1, 1., x2);
  TPad* pad3 = new TPad("pad3", "pad3", 0., 0., 1., x1);
  pad1->SetRightMargin (0.03); pad2->SetRightMargin(pad1->GetRightMargin()); pad3->SetRightMargin(pad1->GetRightMargin());
  pad1->SetLeftMargin  (0.14); pad2->SetLeftMargin (pad1->GetLeftMargin ()); pad3->SetLeftMargin (pad1->GetLeftMargin ());
  pad1->SetBottomMargin(0.03); pad2->SetBottomMargin(0.11); pad3->SetBottomMargin(0.40);
  pad1->SetTopMargin   (0.10); pad2->SetTopMargin   (0.07); pad3->SetTopMargin   (0.03);
  pad1->Draw(); pad2->Draw(); pad3->Draw(); pad2->Draw(); //re-draw pad 2 to put label on white space in pad 2

  // Draw the data and fit components
  pad1->cd();

  //Configure the data style
  gdata->SetMarkerStyle(20);
  gdata->SetMarkerSize(1.3);
  gdata->SetLineWidth(2);

  //Configure the total fit (S+B) style
  htotal->SetLineColor(kRed);
  htotal->SetMarkerColor(kRed);
  htotal->SetFillColor(0);
  htotal->SetFillStyle(3003);
  // htotal->SetMarkerStyle(20);
  htotal->SetMarkerSize(0.);
  htotal->SetLineWidth(3);
  htotal->SetLineColor(kRed);
  htotal->SetTitle("");
  htotal->SetXTitle("");
  if(s_over_sb) htotal->SetYTitle(Form("S/(S+B) weighted Events / %.1f GeV", htotal->GetBinWidth(1)));
  else          htotal->SetYTitle(Form("Events / %.1f GeV", htotal->GetBinWidth(1)));
  htotal->GetYaxis()->SetTitleSize(0.08);
  htotal->GetYaxis()->SetTitleOffset(0.85);
  htotal->GetYaxis()->SetLabelSize(0.065);
  htotal->GetXaxis()->SetLabelSize(0.);

  //Configure the background component style
  hbkg->SetLineColor(kRed);
  hbkg->SetMarkerColor(kRed);
  hbkg->SetLineWidth(3);
  hbkg->SetLineStyle(kDashed);

  //Configure the Z->mumu component style
  if(hzmumu) {
    hzmumu->SetLineColor(kGreen);
    hzmumu->SetMarkerColor(kGreen);
    hzmumu->SetLineWidth(3);
    hzmumu->SetLineStyle(kDashed);
  }

  //Configure the signal component style
  hsignal->SetLineColor(kBlue);
  hsignal->SetMarkerColor(kBlue);
  hsignal->SetLineWidth(3);
  hsignal->SetLineStyle(kDashed);

  //Draw the results
  htotal->Draw("L");
  if(unblind_) {
    hbkg->Draw("L same");
    if(hzmumu) hzmumu->Draw("L same");
    hsignal->Draw("L same");
  }
  gdata->Draw("PZ");
  htotal->GetXaxis()->SetRangeUser(xmin, xmax);


  //Add a legend
  TLegend leg(0.6, 0.5, 0.85, 0.85);
  leg.AddEntry(gdata, "Data", "PE");
  leg.AddEntry(htotal, "Background+signal", "LF");
  if(unblind_) {
    if(zprime) leg.AddEntry(hbkg, "Background", "L");
    else       leg.AddEntry(hbkg, "Parametric background", "L");
    if(hzmumu) leg.AddEntry(hzmumu, "Z#rightarrow#mu#mu", "L");
    if(zprime) leg.AddEntry(hsignal, "Z'#rightarrowe#mu", "L");
    else       leg.AddEntry(hsignal, "Z#rightarrowe#mu", "L");
  }
  leg.SetTextSize(0.065);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetLineColor(0);
  leg.SetLineStyle(0);
  leg.Draw();

  //Make the background subtracted plots
  pad2->cd();

  //Make the data / total fit and data / background component distributions
  TGraphAsymmErrors* gDiff_s = (TGraphAsymmErrors*) gdata->Clone("gDiff_s"); //data / S+B
  TGraphAsymmErrors* gDiff_b = (TGraphAsymmErrors*) gdata->Clone("gDiff_b"); //data / B-only
  TH1* hPull_s = (TH1*) htotal->Clone("hPull_s"); hPull_s->Reset(); //(data - model)/unc (S+B-only)
  TH1* hPull_b = (TH1*) htotal->Clone("hPull_b"); hPull_b->Reset(); //(data - model)/unc (B-only)
  double min_diff(1.e10), max_diff(-1.e10);
  //Debug printout if needed
  if(debug_) {
    printf("Bin:     data        B +-  sigma_B       S+B +- sigma_S+B  pull_B  pull_S+B\n");

  }
  for(int bin = 0; bin < nbins; ++bin) {
    //Retrieve the data point and corresponding model value
    double x, y;
    gdata->GetPoint(bin, x, y);
    const double err_high = gdata->GetErrorYhigh(bin);
    const double err_low  = gdata->GetErrorYlow (bin);
    const double tot_v    = htotal->GetBinContent(bin+1);
    const double tot_e    = htotal->GetBinError  (bin+1);
    const double bkg_v    = hbackground->GetBinContent(bin+1);
    const double bkg_e    = hbackground->GetBinError  (bin+1);

    //Calculate the differences and errors (data-only)
    const double val_s      = (tot_v > 0.) ? y - tot_v : 0.;
    const double val_s_high = (tot_v > 0.) ? err_high  : 0.;
    const double val_s_low  = (tot_v > 0.) ? err_low   : 0.;
    const double val_b      = (bkg_v > 0.) ? y - bkg_v : 0.;
    const double val_b_high = (bkg_v > 0.) ? err_high  : 0.;
    const double val_b_low  = (bkg_v > 0.) ? err_low   : 0.;

    //Calculate the pulls, combining the data and model errors
    const double data_err_s = (y > tot_v) ? err_low : err_high;
    const double data_err_b = (y > bkg_v) ? err_low : err_high;
    // const double err_s      = (tot_e <= 0.) ? 0. : (err_mode_ == 0) ? sqrt(tot_e*tot_e + data_err_s*data_err_s) : sqrt(data_err_s*data_err_s - tot_e*tot_e);
    // const double err_b      = (bkg_e <= 0.) ? 0. : (err_mode_ == 0) ? sqrt(bkg_e*bkg_e + data_err_b*data_err_b) : sqrt(data_err_b*data_err_b - bkg_e*bkg_e);
    // if(!std::isfinite(err_s) || !std::isfinite(err_b)) printf(" %s: %s bin %i has non-finite error(s): err_s = %f; err_b = %f\n",
    //                                                           __func__, tag.Data(), bin, err_s, err_b);
    // const double pull_s     = (tot_e > 0.) ? (y - tot_v) / err_s : 0.;
    // const double pull_b     = (bkg_e > 0.) ? (y - bkg_v) / err_b : 0.;

    const double pull_s     = (tot_v > 0.) ? (y - tot_v) / sqrt(tot_v) : 0.;
    const double pull_b     = (bkg_v > 0.) ? (y - bkg_v) / sqrt(bkg_v) : 0.;

    //Set the points
    gDiff_s->SetPoint      (bin, x, val_s);
    gDiff_s->SetPointEYhigh(bin, val_s_high);
    gDiff_s->SetPointEYlow (bin, val_s_low );
    gDiff_b->SetPoint      (bin, x, val_b);
    gDiff_b->SetPointEYhigh(bin, val_b_high);
    gDiff_b->SetPointEYlow (bin, val_b_low );
    hPull_s->SetBinContent(bin+1, pull_s);
    hPull_b->SetBinContent(bin+1, pull_b);

    max_diff = max(max_diff, val_b);
    min_diff = min(min_diff, val_b);

    //Debug printout if needed
    if(debug_) {
      printf(" %2i: %8.0f %8.1f +- %8.2f  %8.1f +- %8.2f   %5.2f    %5.2f\n",
             bin+1, y, bkg_v, bkg_e, tot_v, tot_e, pull_b, pull_s);
    }
  }

  auto gDiff = (unblind_) ? gDiff_b : gDiff_s; //Difference plot to use, data-total if blinded, data-bkg-only otherwise

  //Make a background uncertainty plot
  TH1* hBkg_unc = (TH1*) hbackground->Clone("hBkg_unc");
  for(int ibin = 1; ibin <= hbackground->GetNbinsX(); ++ibin) {
    if(hbackground->GetBinContent(ibin) <= 0.01) hBkg_unc->SetBinContent(ibin, 0.);
    else {
      hBkg_unc->SetBinContent(ibin, 1.);
      if(hbackground->GetBinError(ibin) > 0.) {
        hBkg_unc->SetBinError(ibin, hbackground->GetBinError(ibin) / hbackground->GetBinContent(ibin));
      } else {
        hBkg_unc->SetBinError(ibin, 0.1 / hbackground->GetBinContent(ibin));
      }
    }
  }
  hBkg_unc->SetFillColor(kGray+3);
  hBkg_unc->SetLineColor(kGray+3);
  hBkg_unc->SetFillStyle(3003);

  //Add the signal fit to the difference plot
  hzemu->SetFillColor(0);
  hzemu->SetLineWidth(3);
  hzemu->SetLineStyle(kSolid);

  //Draw the results
  const double rmin(min_diff - 0.20*(max_diff-min_diff)), rmax(max_diff + 0.20*(max_diff-min_diff));
  hBkg_unc->Draw("E2");
  if(unblind_) {
    hzemu->Draw("L same");
    hzemu->GetYaxis()->SetRangeUser(rmin, rmax); //necessary due to y-axis importing from clone
  }
  gDiff->Draw("PZ");
  hBkg_unc->GetXaxis()->SetRangeUser(xmin, xmax);

  //Add a reference line for perfect agreement
  TLine* line = new TLine(xmin, 0., xmax, 0.);
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  //Configure the titles and axes
  float txt_scale = (1.-x2)/(x2-x1);
  hBkg_unc->GetYaxis()->SetRangeUser(rmin, rmax);
  hBkg_unc->GetYaxis()->SetNdivisions(505);
  hBkg_unc->SetTitle("");
  hBkg_unc->SetXTitle("");
  hBkg_unc->GetXaxis()->SetLabelSize(0.);
  hBkg_unc->GetYaxis()->SetLabelSize(txt_scale*htotal->GetYaxis()->GetLabelSize());
  hBkg_unc->GetYaxis()->SetTitleSize(txt_scale*htotal->GetYaxis()->GetTitleSize());
  hBkg_unc->GetYaxis()->SetTitleOffset(0.22);
  if(unblind_) hBkg_unc->SetYTitle("Data-Bkg");
  else         hBkg_unc->SetYTitle("Data-Fit");


  //Make a pull plot
  pad3->cd();

  txt_scale = (1.-x2)/x1;
  auto hPull = hPull_s; //FIXME: Decide whether pulls should use background or total model
  hPull->SetLineColor(kAtlantic);
  hPull->SetFillColor(kAtlantic);
  hPull->SetFillStyle(1000);
  hPull->Draw("hist");
  hPull->GetYaxis()->SetRangeUser(-3,3);
  hPull->GetYaxis()->SetNdivisions(505);
  hPull->SetTitle("");
  hPull->SetXTitle("m_{e#mu} [GeV]");
  hPull->SetYTitle("#frac{Data-Fit}{#sigma_{Fit}}");
  hPull->GetYaxis()->SetLabelSize(txt_scale*htotal->GetYaxis()->GetLabelSize());
  hPull->GetYaxis()->SetTitleSize(txt_scale*htotal->GetYaxis()->GetTitleSize());
  hPull->GetXaxis()->SetLabelSize(hPull->GetYaxis()->GetLabelSize());
  hPull->GetXaxis()->SetTitleSize(hPull->GetYaxis()->GetTitleSize());
  hPull->GetXaxis()->SetLabelOffset(0.01);
  hPull->GetYaxis()->SetLabelOffset(0.015);
  hPull->GetXaxis()->SetTitleOffset(0.81);
  hPull->GetYaxis()->SetTitleOffset(0.27);

  //Add a reference line for perfect agreement
  TLine* line_2 = new TLine(xmin, 0., xmax, 0.);
  line_2->SetLineColor(kBlack);
  line_2->SetLineWidth(2);
  line_2->SetLineStyle(kDashed);
  line_2->Draw("same");

  //Add the CMS label
  pad1->cd();
  draw_cms_label(pad1->GetLeftMargin());
  draw_luminosity();
  draw_category(tag, zprime, pad1->GetLeftMargin());

  //Print a linear and a log version of the distribution
  const double min_signal = hmin(hzemu, -1.e10);
  cout << "min signal = " << min_signal << endl;
  const double min_val = std::max(0.1, std::min(gmin(gdata), hmin(htotal)));
  const double max_val = std::max(gmax(gdata), hmax(htotal));
  htotal->GetYaxis()->SetRangeUser((min_signal < 0.) ? 1.1*min_signal : 0., 1.3*(max_val + 1.05*sqrt(max_val)));
  c->SaveAs(Form("%s%s.%s", outdir.Data(), tag.Data(), file_type_.Data()));
  c->SaveAs(Form("%s%s.root", outdir.Data(), tag.Data()));
  const double plot_min = std::min(std::max(0.2, 0.2*hmax(hsignal)), 0.2*min_val);
  const double plot_max = plot_min*pow(10, 1.7*log10(max_val/plot_min));
  htotal->GetYaxis()->SetRangeUser(plot_min, plot_max);
  pad1->SetLogy();
  c->SaveAs(Form("%s%s_logy.%s", outdir.Data(), tag.Data(), file_type_.Data()));
  c->SaveAs(Form("%s%s_logy.root", outdir.Data(), tag.Data()));

  //Clean up after printing
  delete c;
  delete gDiff_s;
  delete gDiff_b;
  delete hPull_s;
  delete hPull_b;
  delete hBkg_unc;
  delete line;
  delete line_2;

  return 0;
}

//------------------------------------------------------------------------------------------
// Print the fit results for each category in a fit configuration directory
int print_dir(TDirectoryFile* dir, TString tag, TString outdir) {
  int status(0);
  if(!dir) return 1;

  //List of categories
  TList* list = dir->GetListOfKeys();
  if(!list) return 10;
  vector<TDirectoryFile*> subdirs;
  bool subdir(false); //whether there are sub-directories or not
  for(TObject* o : *list) {
    TObject* obj = dir->Get(o->GetName());
    if(!obj) continue;

    //Check if this object is a directory with categories, if so recursively process it
    bool isdir = obj->InheritsFrom(TDirectoryFile::Class());
    if(isdir) {
      auto next_dir = (TDirectoryFile*) obj;
      status += print_dir(next_dir, (tag + "_") + obj->GetName(), outdir);
      subdir = true;
      if(do_single_) return status;
      subdirs.push_back(next_dir);
    }
  }

  //If this directory doesn't contain a sub-directory, print the histograms within the category
  if(!subdir) { //histogram directory
    if(!only_s_fit_ || tag.Contains("fit_s"))
      status += print_hist({dir}, tag, outdir);
  } else if(do_sig_wt_) { //print the S/(S+B) weighted histogram
    tag += "_s_sb_wt";
    if(is_prelim_) tag += "_prelim";
    status += print_hist(subdirs, tag, outdir, true);
  }
  return status;
}

//------------------------------------------------------------------------------------------
// Print all fit figures
int plot_bemu(TString fname, TString outdir = "figures", bool unblind = false) {
  unblind_ = unblind;

  //Get the fit file
  TFile* file = TFile::Open(fname.Data(), "READ");
  if(!file) return 1;

  TDirectoryFile *prefit = (TDirectoryFile*) file->Get("shapes_prefit");
  TDirectoryFile *fit_b  = (TDirectoryFile*) file->Get("shapes_fit_b");
  TDirectoryFile *fit_s  = (TDirectoryFile*) file->Get("shapes_fit_s");

  if(!prefit || (unblind && (!fit_b || !fit_s))) {
    cout << "Fit directories not found!\n";
    return 2;
  }

  //Create the figure directory
  if(!outdir.EndsWith("/")) outdir += "/";
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", outdir.Data(), outdir.Data()));

  //Print the fit configurations: Pre-fit, background-only fit, and background+signal fit
  int status(0);
  if(!do_single_) status += print_dir(prefit, "prefit", outdir);
  if(!do_single_ && fit_b) status += print_dir(fit_b , "fit_b" , outdir);
  if(fit_s) status += print_dir(fit_s , "fit_s" , outdir);

  cout << "Plotting status = " << status << endl;
  return status;
}
