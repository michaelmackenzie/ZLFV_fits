//Draw workspace PDFs
const bool blind_data_ = true;

void draw_cms_label(const bool inside = false, const bool is_prelim = true) {
  TText cmslabel;
  cmslabel.SetNDC();
  cmslabel.SetTextColor(1);
  cmslabel.SetTextSize(0.06);
  cmslabel.SetTextAlign(11);
  cmslabel.SetTextAngle(0);
  cmslabel.SetTextFont(61);
  const float label_y(0.915f);
  if(inside) cmslabel.DrawText(0.16, 0.82   , "CMS");
  else       cmslabel.DrawText(0.10, label_y, "CMS");
  if(is_prelim) {
    cmslabel.SetTextFont(52);
    cmslabel.SetTextSize(0.76*cmslabel.GetTextSize());
    if(inside) cmslabel.DrawText(0.16, 0.77   , "Preliminary");
    else       cmslabel.DrawText(0.20, label_y, "Preliminary");
  }
}

void draw_lumi_label() {
  TLatex lumilabel;
  lumilabel.SetNDC();
  lumilabel.SetTextFont(42);
  lumilabel.SetTextSize(0.045);
  lumilabel.SetTextAlign(31);
  lumilabel.SetTextAngle(0);
  const float label_y(0.915f);
  const int year = -1; //default to Run 2
  TString period = (year > 2000) ? Form("%i, ", year) : "";
  const double lum = (year == 2016) ? 36.33 : (year == 2017) ? 41.48 : (year == 2018) ? 59.83 : 137.64;
  lumilabel.DrawLatex(0.97, label_y, Form("%s%.0f fb^{-1} (13 TeV)",period.Data(),lum));
}

void draw_bdt_region(TString bin) {
  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.045);
  label.SetTextAlign(11);
  label.SetTextAngle(0);
  float bdt_low(0.f), bdt_high(1.f);
  if     (bin == "bin1") {bdt_low = 0.3f; bdt_high = 0.7f;}
  else if(bin == "bin2") {bdt_low = 0.7f; bdt_high = 0.9f;}
  else if(bin == "bin3") {bdt_low = 0.9f; bdt_high = 1.0f;}
  label.DrawLatex(0.16, 0.77, Form("%.1f < BDT < %.1f", bdt_low, bdt_high));
}

TString get_pdf_title(TString pdf_name) {
  TString title = "";
  if(pdf_name.Contains("_gs"))   title += ""; //"Gaussian + ";
  if(pdf_name.Contains("pol1"))  title += "Poly (1)";
  if(pdf_name.Contains("pol2"))  title += "Poly (2)";
  if(pdf_name.Contains("pol3"))  title += "Poly (3)";
  if(pdf_name.Contains("pol4"))  title += "Poly (4)";
  if(pdf_name.Contains("cheb1")) title += "Chebychev (1)";
  if(pdf_name.Contains("cheb2")) title += "Chebychev (2)";
  if(pdf_name.Contains("cheb3")) title += "Chebychev (3)";
  if(pdf_name.Contains("cheb4")) title += "Chebychev (4)";
  if(pdf_name.Contains("cheb5")) title += "Chebychev (5)";
  if(pdf_name.Contains("plaw1")) title += "Power (1)";
  if(pdf_name.Contains("plaw2")) title += "Power (2)";
  if(pdf_name.Contains("plaw3")) title += "Power (3)";
  if(pdf_name.Contains("exp1"))  title += "Expo (1)";
  if(pdf_name.Contains("exp2"))  title += "Expo (2)";
  if(pdf_name.Contains("exp3"))  title += "Expo (3)";
  return title;
}

void set_axis_styles(TAxis* upper_x, TAxis* upper_y,
                     TAxis* lower_x, TAxis* lower_y) {

  TGaxis::SetMaxDigits(3);

  upper_x->SetTitle("");
  upper_y->SetTitle(Form("Events / %.1f GeV", upper_x->GetBinWidth(1)));
  lower_y->SetTitle("#frac{N_{data} - N_{fit}}{#sigma_{data}}");
  lower_x->SetTitle("m_{e#mu} [GeV]");

  upper_x->SetLabelSize(0.);
  upper_y->SetLabelSize(0.05);
  upper_y->SetTitleSize(0.06);
  upper_y->SetTitleOffset(1.07);

  lower_y->SetRangeUser(-3.99, 3.99);
  lower_y->SetNdivisions(205);
  lower_y->SetLabelSize(0.11);
  lower_y->SetLabelOffset(0.01);
  lower_y->SetTitleSize(0.13);
  lower_y->SetTitleOffset(0.41);
  lower_x->SetLabelSize(0.11);
  lower_x->SetLabelOffset(0.008);
  lower_x->SetTitleSize(0.13);
  lower_x->SetTitleOffset(0.90);

  upper_x->SetRangeUser(70., 110.);
  lower_x->SetRangeUser(70., 110.);
}

TGraphAsymmErrors* envelope_band(RooMultiPdf* multipdf, TH1* hdata, TH1* hzmm, RooRealVar* obs,
                                 TGraphAsymmErrors *& data_errors, TGraphAsymmErrors *& data_sig_errors,
                                 TGraphAsymmErrors *& band_errors) {
  vector<TH1*> hpdfs;
  const int npdfs = multipdf->getNumPdfs();
  const float blind_min(86.f), blind_max(96.f);
  if(npdfs <= 0) return nullptr;
  for(int ipdf = 0; ipdf < npdfs; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    hpdfs.push_back(pdf->createHistogram(Form("pdf_hist_%i", ipdf), *obs));
  }
  const double scale = (hdata->Integral() - hzmm->Integral()) / hpdfs[0]->Integral();
  const int nbins = hpdfs[0]->GetNbinsX();
  double xvals[nbins], yvals[nbins], xerrs[nbins], yerr_ups[nbins], yerr_downs[nbins]; //band
  double pulls[nbins], pulls_sig[nbins], pull_yerrs[nbins]; //data pulls
  double rvals[nbins], ryerr_ups[nbins], ryerr_downs[nbins]; //band ratio
  for(int index = 0; index < nbins; ++index) {
    double ymin(1.e10), ymax(-1.e10), x(hpdfs[0]->GetBinCenter(index+1)), ndata(hdata->GetBinContent(index+1));
    for(int ipdf = 0; ipdf < npdfs; ++ipdf) {
      const double y = hpdfs[ipdf]->GetBinContent(index+1)*scale + hzmm->GetBinContent(index+1);
      ymin = min(ymin, y);
      ymax = max(ymax, y);
    }
    bool blind = blind_data_ && x >= blind_min && x <= blind_max;
    const double y = (ymax + ymin)/2.;
    xvals[index] = x;
    yvals[index] = y;
    yerr_ups[index] = (ymax - ymin)/2.;
    yerr_downs[index] = (ymax - ymin)/2.;
    xerrs[index] = 0.; //no bin widths
    if(blind) {
      pulls_sig[index] = (ndata - y) / sqrt(ndata);
      pulls    [index] = -999.;
    } else {
      pulls    [index] = (ndata - y) / sqrt(ndata);
      pulls_sig[index] = -999.;
    }
    pull_yerrs[index] = 1.;
    rvals[index] = 0.;
    ryerr_ups  [index] = (y > 0.) ? (yerr_ups  [index]) / sqrt(y) : 0.;
    ryerr_downs[index] = (y > 0.) ? (yerr_downs[index]) / sqrt(y) : 0.;
  }
  data_errors     = new TGraphAsymmErrors(nbins, xvals, pulls    , xerrs, xerrs, pull_yerrs, pull_yerrs);
  data_sig_errors = new TGraphAsymmErrors(nbins, xvals, pulls_sig, xerrs, xerrs, pull_yerrs, pull_yerrs);
  band_errors     = new TGraphAsymmErrors(nbins, xvals, rvals    , xerrs, xerrs, ryerr_ups, ryerr_downs);
  return new TGraphAsymmErrors(nbins, xvals, yvals, xerrs, xerrs, yerr_ups, yerr_downs);
}

int debug_zemu_ws(const char* bin = "bin3") {

  // Get the name of the observable
  TString bin_s(bin);
  TString obs_name = "lepm_1";
  if(bin_s.BeginsWith("bin")) obs_name = obs_name + bin_s[bin_s.Sizeof()-2];
  else obs_name = "lepm_" + bin_s;

  ////////////////////////////////////////
  // Retrieve the input background PDFs
  ////////////////////////////////////////

  TFile* f_bkg = TFile::Open(Form("workspace_mk2bkg_v1_%s.root", bin), "READ");
  if(!f_bkg) return 1;
  auto ws = (RooWorkspace*) f_bkg->Get("workspace_background");
  if(!ws) return 2;
  auto multipdf = (RooMultiPdf*) ws->pdf(Form("multipdf_%s", bin));
  if(!multipdf) return 3;
  auto cat = (RooCategory*) ws->obj(Form("pdfindex_%s", bin));
  if(!cat) return 4;
  cout << "Using obs name " << obs_name.Data() << endl;
  auto obs = (RooRealVar*) ws->var(obs_name);
  if(!obs) return 5;
  const float xmin(70.f), xmax(110.f), blind_min(86.f), blind_max(96.f);
  obs->setRange("full", xmin, xmax);
  obs->setRange("LowSideband", xmin, blind_min-1.e-4);
  obs->setRange("HighSideband", blind_max+1.e-4, xmax);
  obs->setRange("BlindRegion", blind_min, blind_max);
  obs->setUnit("GeV");

  ////////////////////////////////////////
  // Retrieve the input data
  ////////////////////////////////////////

  // workspace_mk2dataAndMC_v1_bin1.root workspace_data:real_data_bin1
  TFile* f_data = TFile::Open(Form("workspace_mk2dataAndMC_v1_%s.root", bin), "READ");
  if(!f_data) return 11;
  auto ws_data = (RooWorkspace*) f_data->Get("workspace_data");
  if(!ws_data) return 12;
  auto data = ws_data->data(Form("real_data_%s", bin));
  if(!data) return 13;

  ////////////////////////////////////////
  // Retrieve the Z->mumu PDF
  ////////////////////////////////////////

  // workspace_mk2zmm_v1_bin1.root workspace_zmm:zmm_dcb_pdf_bin1
  TFile* f_zmm = TFile::Open(Form("workspace_mk2zmm_v1_%s.root", bin), "READ");
  if(!f_zmm) return 21;
  auto ws_zmm = (RooWorkspace*) f_zmm->Get("workspace_zmm");
  if(!ws_zmm) return 22;
  auto zmm = ws_zmm->pdf(Form("zmm_dcb_pdf_%s", bin));
  if(!zmm) return 23;
  auto zmm_norm = (RooRealVar*) ws_zmm->var(Form("zmm_dcb_pdf_%s_norm", bin));
  if(!zmm_norm) return 24;

  ////////////////////////////////////////
  // Make a plot of the inputs
  ////////////////////////////////////////

  //Create a legend
  TLegend* leg = new TLegend(0.11, 0.65, 0.94, 0.89);
  leg->SetNColumns(3);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  auto frame = obs->frame();
  data->plotOn(frame, RooFit::Name("data"));

  const int ncats = cat->numTypes();
  // const int colors[] = {kRed-6, kGreen-6, kOrange-4, kViolet-2};
  const int styles[] = {kSolid, 6, 10};
  const int colors[] = {kRed-6, kGreen-6, kOrange-4, kViolet-2, kYellow-3, kAtlantic, kOrange+1, kSpring+5};
  const int ncolors = sizeof(colors)/sizeof(*colors);
  const int nstyles = sizeof(styles)/sizeof(*styles);
  vector<TString> titles;
  for(int ipdf = 0; ipdf < ncats; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    // if(ipdf % 2 == 0) continue;
    TString title = get_pdf_title(pdf->GetName());
    titles.push_back(title);
    pdf->plotOn(frame, RooFit::LineColor(colors[ipdf % ncolors]), RooFit::LineStyle(ipdf % 10 + 3), //styles[min(nstyles-1, ipdf / 2)]),
                RooFit::Title(title), RooFit::Name(Form("pdf_%i", ipdf)));
  }

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  TCanvas c; c.SetRightMargin(0.05);
  frame->Draw();
  leg->AddEntry("data", "Data", "PL");
  for(int ipdf = 0; ipdf < ncats; ++ipdf) leg->AddEntry(Form("pdf_%i", ipdf), titles[ipdf], "L");
  leg->Draw("same");
  frame->GetYaxis()->SetRangeUser(0.1, 1.4*frame->GetMaximum());
  c.SaveAs(Form("inputs_%s.pdf", bin));
  c.SaveAs(Form("inputs_%s.root", bin));

  // Draw the Z->mumu PDF
  frame = obs->frame();
  zmm->plotOn(frame, RooFit::Normalization(zmm_norm->getVal()));
  frame->Draw();
  c.SaveAs(Form("inputs_zmm_%s.pdf", bin));
  c.SaveAs(Form("inputs_zmm_%s.root", bin));

  ////////////////////////////////////////
  // Re-fit the Z->mumu+PDFs to the data
  ////////////////////////////////////////

  // Blind the data for this plot
  TH1* h_data     = data->createHistogram("hdata"    , *obs);
  TH1* h_sig_data = data->createHistogram("hdata_sig", *obs); //data in the signal region
  if(blind_data_) {
    const int blind_low_bin  = h_data->FindBin(blind_min+0.01);
    const int blind_high_bin = h_data->FindBin(blind_max-0.01);
    for(int ibin = 1; ibin <= h_data->GetNbinsX(); ++ibin) {
      if(ibin >= blind_low_bin && ibin <= blind_high_bin) {
        h_data->SetBinContent(ibin, 0.);
        h_data->SetBinError(ibin, 0.);
      } else {
        h_sig_data->SetBinContent(ibin, 0.);
        h_sig_data->SetBinError(ibin, 0.);
      }
    }
  }
  auto blind_data  = new RooDataHist("blind_data", "Blinded data", *obs, h_data);
  auto signal_data = new RooDataHist("signal_data", "Signal region data", *obs, h_sig_data);

  TCanvas c2("c2","c2", 700, 800);
  TPad pad1("pad1","pad1",0.,0.3,1.,1.0); pad1.Draw();
  TPad pad2("pad2","pad2",0.,0.0,1.,0.3); pad2.Draw();
  pad1.SetTopMargin(0.1);
  pad1.SetRightMargin(0.03);
  pad1.SetLeftMargin(0.132);
  pad1.SetBottomMargin(0.02);
  pad2.SetTopMargin(0.03);
  pad2.SetRightMargin(pad1.GetRightMargin());
  pad2.SetLeftMargin(pad1.GetLeftMargin());
  pad2.SetBottomMargin(0.26);

  pad1.cd();
  frame = obs->frame();
  auto frame2 = obs->frame();
  blind_data->plotOn(frame, RooFit::Name("data"));
  signal_data->plotOn(frame, RooFit::Name("signal_data"), RooFit::MarkerColor(kBlue), RooFit::LineColor(kBlue));
  const double frac_zmm = zmm_norm->getVal() / data->sumEntries();
  RooRealVar* pdf_frac = new RooRealVar("pdf_frac", "PDF fraction", 1. - frac_zmm);
  RooRealVar* zmm_val = new RooRealVar("zmm_val", "Z->mumu events", zmm_norm->getVal());
  pdf_frac->setConstant(true);
  zmm_val->setConstant(true);
  for(int ipdf = 0; ipdf < ncats; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    // RooAddPdf* tot_pdf = new RooAddPdf(Form("tot_pdf_%i", ipdf), titles[ipdf], *pdf, *zmm, *pdf_frac);
    RooRealVar* pdf_val = new RooRealVar(Form("pdf_val_%i", ipdf), "PDF events", data->sumEntries() - zmm_norm->getVal());
    RooAddPdf* tot_pdf = new RooAddPdf(Form("tot_pdf_%i", ipdf), titles[ipdf], RooArgList(*pdf, *zmm), RooArgList(*pdf_val, *zmm_val));
    tot_pdf->fitTo(*data, RooFit::PrintLevel(-1), RooFit::Warnings(0), RooFit::PrintEvalErrors(-1),
                   RooFit::Range("LowSideband,HighSideband"), RooFit::NormRange("LowSideband,HighSideband"));
    tot_pdf->plotOn(frame, RooFit::LineColor(colors[ipdf % ncolors]), RooFit::LineStyle(ipdf % 10 + 3),
                    //RooFit::LineStyle((ipdf >= ncolors) ? kDashed : kSolid),
                    RooFit::Title(titles[ipdf]), RooFit::Name(Form("pdf_%i", ipdf)), RooFit::Range("full"), RooFit::NormRange("LowSideband,HighSideband"));
    if(ipdf == 0) {
      tot_pdf->plotOn(frame, RooFit::Components(zmm->GetName()), RooFit::LineColor(kGreen), RooFit::Name("zmm"),
                      RooFit::Range("full"), RooFit::NormRange("LowSideband,HighSideband"));
    }
    auto h = frame->pullHist("data", Form("pdf_%i", ipdf));
    h->SetLineWidth(2);
    h->SetLineColor(colors[ipdf % ncolors]);
    h->SetLineStyle(ipdf % 10 + 3);
    h->SetMarkerColor(h->GetLineColor());
    h->SetMarkerStyle(2);
    h->SetMarkerSize(1.0);
    frame2->addPlotable(h, "PE1");
  }
  frame->Draw();
  leg->AddEntry("zmm", "Z->#mu#mu", "L");
  leg->SetTextSize(0.055);
  leg->Draw("same");
  frame->GetYaxis()->SetRangeUser(0.1, 1.4*frame->GetMaximum());

  pad2.cd();
  frame2->Draw();

  frame->SetTitle("");
  frame2->SetTitle("");

  // TGaxis::SetExponentOffset(-0.05, 0.01, "Y");
  set_axis_styles(frame->GetXaxis(), frame->GetYaxis(),
                  frame2->GetXaxis(), frame2->GetYaxis());

  TLine line(frame->GetXaxis()->GetXmin(), 0., frame->GetXaxis()->GetXmax(), 0.);
  line.SetLineWidth(2);
  line.SetLineColor(kBlack);
  line.SetLineStyle(kDashed);
  line.Draw("same");
  pad2.SetGridy();

  //CMS prelim drawing
  pad1.cd();
  draw_cms_label(true, false);
  draw_lumi_label();

  c2.SaveAs(Form("inputs_refit_%s.pdf", bin));
  c2.SaveAs(Form("inputs_refit_%s.root", bin));

  // Create a envelope band from the input functions
  auto h_zmm = zmm->createHistogram("h_zmm", *obs);
  h_zmm->Scale(zmm_norm->getVal() / h_zmm->Integral());
  TGraphAsymmErrors *data_errors(nullptr), *data_sig_errors(nullptr), *band_errors(nullptr);
  auto band = envelope_band(multipdf, data->createHistogram("hdata_2", *obs), h_zmm, obs, data_errors, data_sig_errors, band_errors);
  band->SetLineWidth(1);
  band->SetLineColor(kRed);
  band->SetMarkerSize(0.);
  band->SetFillColor(kGray);

  for(int ibin = 1; ibin <= h_zmm->GetNbinsX(); ++ibin) h_zmm->SetBinError(ibin, 0.);
  h_zmm->SetLineColor(kGreen+2);
  // h_zmm->SetLineStyle(kDashed);
  h_zmm->SetLineWidth(2);
  h_zmm->SetFillColor(0);
  h_zmm->SetFillStyle(0);

  pad1.cd();
  h_data->Draw("EX0");
  h_sig_data->Draw("EX0 SAME");
  band->Draw("E3");
  band->Draw("LX");
  h_zmm->Draw("L0 same");

  h_data->SetLineWidth(2);
  h_data->SetLineColor(kBlack);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.8);
  h_data->SetTitle("");
  h_data->GetXaxis()->SetLabelSize(0.);
  h_data->SetYTitle(Form("Events / %.1f GeV", h_data->GetBinWidth(1)));
  h_data->GetYaxis()->SetRangeUser(0., 1.3*h_data->GetMaximum());

  h_sig_data->SetLineWidth  (h_data->GetLineWidth  ());
  h_sig_data->SetMarkerStyle(h_data->GetMarkerStyle());
  h_sig_data->SetMarkerSize (h_data->GetMarkerSize ());
  h_sig_data->SetLineColor  (kBlue);
  h_sig_data->SetMarkerColor(kBlue);

  // Ratio plot
  pad2.cd();

  if(!data_errors || !band_errors) return 10;
  data_errors->SetTitle("");
  data_errors->SetLineWidth  (h_data->GetLineWidth  ());
  data_errors->SetLineColor  (h_data->GetLineColor  ());
  data_errors->SetMarkerStyle(h_data->GetMarkerStyle());
  data_errors->SetMarkerSize (h_data->GetMarkerSize ());
  data_errors->SetMarkerColor(h_data->GetMarkerColor());
  data_errors->Draw("APEZ");

  data_sig_errors->SetTitle("");
  data_sig_errors->SetLineWidth  (h_sig_data->GetLineWidth  ());
  data_sig_errors->SetLineColor  (h_sig_data->GetLineColor  ());
  data_sig_errors->SetMarkerStyle(h_sig_data->GetMarkerStyle());
  data_sig_errors->SetMarkerSize (h_sig_data->GetMarkerSize ());
  data_sig_errors->SetMarkerColor(h_sig_data->GetMarkerColor());

  band_errors->SetMarkerSize(0.);
  band_errors->SetFillColor(kGray);
  band_errors->Draw("E3");
  data_errors->Draw("PEZ");
  data_sig_errors->Draw("PEZ");

  TGaxis::SetExponentOffset(-0.05, 0.01, "Y");
  set_axis_styles(h_data->GetXaxis(), h_data->GetYaxis(),
                  data_errors->GetXaxis(), data_errors->GetYaxis());

  //CMS prelim drawing
  pad1.cd();
  draw_cms_label(true, false);
  draw_lumi_label();
  draw_bdt_region(bin);

  leg = new TLegend(0.60, 0.55, 0.94, 0.89);
  leg->SetNColumns(1);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.055);

  leg->AddEntry(h_data, "Data", "PE");
  leg->AddEntry(band, "Background", "LF");
  // leg->AddEntry(band_errors, "Envelope", "F");
  leg->AddEntry(h_zmm, "Z#rightarrow#mu#mu", "L");
  leg->Draw();

  c2.SaveAs(Form("inputs_refit_band_%s.pdf", bin));
  c2.SaveAs(Form("inputs_refit_band_%s.root", bin));
  return 0;
}
