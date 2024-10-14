//Draw workspace PDFs

bool blind_ = true;

TString get_pdf_title(TString pdf_name) {
  TString title = "";
  if(pdf_name.Contains("_gs"))   title += "Gaussian + ";
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

int debug_zprime_ws(TString card = "datacard_zprime_bdt_0d3_0d7_v13_mass-110.0_mp0.txt") {

  // Get the name of the observable
  TString obs_name = "mass_ll";
  const char* bin = (card.Contains("0d3_0d7")) ? "bin1" : "bin2";

  ////////////////////////////////////////
  // Retrieve the input background PDFs
  ////////////////////////////////////////

  TObjArray *tx = card.Tokenize("_");
  TString mp = ((TObjString*) (tx->At(tx->GetEntries()-1)))->String();
  mp.ReplaceAll(".txt", "");
  TString mass_s = ((TObjString*) (tx->At(tx->GetEntries()-2)))->String();
  mass_s.ReplaceAll("mass-", "");
  const float mass = std::stof(mass_s.Data());

  TString ws_name = "WorkspaceScanBKG/" + card;
  ws_name.ReplaceAll("datacard_zprime", "workspace_scanbkg_v2");
  ws_name.ReplaceAll("_mass-"+mass_s, "");
  ws_name.ReplaceAll(".txt", ".root");

  TFile* f_bkg = TFile::Open(ws_name, "READ");
  if(!f_bkg) return 1;
  auto ws = (RooWorkspace*) f_bkg->Get("ws_bkg");
  if(!ws) return 2;
  auto multipdf = (RooMultiPdf*) ws->pdf(Form("multipdf_%s", bin));
  if(!multipdf) {ws->Print(); return 3;}
  auto cat = (RooCategory*) ws->obj(Form("pdfindex_%s", bin));
  if(!cat) {ws->Print(); return 4;}
  cout << "Using obs name " << obs_name.Data() << endl;
  auto obs = (RooRealVar*) ws->var(obs_name);
  if(!obs) {ws->Print(); return 5;}
  auto data = (RooDataHist*) ws->data("data_obs");
  if(!data) return 6;
  const float xmin  = obs->getMin();
  const float xmax  = obs->getMax();
  const float sigma = 0.02*mass;
  const int nbin_data = (xmax-xmin)/(0.5*sigma); //0.5 sigma binning
  obs->setBins(nbin_data);

  // Create a blinded data version
  TH1* h_data = data->createHistogram("hdata", *obs);
  const float blind_min = h_data->GetBinLowEdge(h_data->FindBin(mass - sigma + 1.e-4));
  const float blind_max = h_data->GetBinLowEdge(h_data->FindBin(mass + sigma - 1.e-4)) + h_data->GetBinWidth(1);
  if(nbin_data != h_data->GetNbinsX()) {
    cout << "Binning disagreement! N(bins) = " << nbin_data << ", data = " << h_data->GetNbinsX() << endl;
  }
  for(int ibin = h_data->FindBin(blind_min+0.01); ibin <= h_data->FindBin(blind_max-0.1); ++ibin) {
    h_data->SetBinContent(ibin, 0.);
    h_data->SetBinError(ibin, 0.);
  }
  obs->setRange("full", xmin, xmax);
  obs->setRange("LowSideband", xmin, blind_min-1.e-4);
  obs->setRange("HighSideband", blind_max+1.e-4, xmax);
  obs->setRange("BlindRegion", blind_min, blind_max);
  auto blind_data = new RooDataHist("blind_data", "Blinded data", *obs, h_data);
  auto plot_data = (blind_) ? blind_data : data;

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
  plot_data->plotOn(frame, RooFit::Name("data"));
  // if(blind_) data->plotOn(frame, RooFit::Name("unblinded_data"), RooFit::Invisible());

  const int ncats = cat->numTypes();
  const int colors[] = {kRed-6, kGreen-6, kOrange-4, kViolet-4};
  // const int colors[] = {kRed-6, kGreen-6, kOrange-4, kViolet-4, kYellow-3, kAtlantic, kOrange+1, kSpring+5};
  const int ncolors = sizeof(colors)/sizeof(*colors);
  vector<TString> titles;
  for(int ipdf = 0; ipdf < ncats; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    TString title = get_pdf_title(pdf->GetName());
    titles.push_back(title);
    pdf->plotOn(frame, RooFit::LineColor(colors[ipdf % ncolors]), RooFit::LineStyle((ipdf >= ncolors) ? kDashed : kSolid),
                RooFit::Title(title), RooFit::Name(Form("pdf_%i", ipdf)), RooFit::NormRange((blind_) ? "LowSideband,HighSideband" : "full"));
  }

  TCanvas c; c.SetRightMargin(0.05);
  frame->Draw();
  leg->AddEntry("data", "Data", "PL");
  for(int ipdf = 0; ipdf < ncats; ++ipdf) leg->AddEntry(Form("pdf_%i", ipdf), titles[ipdf], "L");
  leg->Draw("same");
  frame->GetYaxis()->SetRangeUser(0.1, 1.4*frame->GetMaximum());
  c.SaveAs(Form("inputs_%s_%s.png", mp.Data(), bin));

  ////////////////////////////////////////
  // Re-fit the PDFs to the data
  ////////////////////////////////////////

  TCanvas c2("c2","c2", 700, 800);
  TPad pad1("pad1","pad1",0.,0.4,1.,1.0); pad1.Draw();
  TPad pad2("pad2","pad2",0.,0.0,1.,0.4); pad2.Draw();
  pad1.SetTopMargin(0.1);
  pad1.SetRightMargin(0.05);
  pad1.SetBottomMargin(0.01);
  pad2.SetTopMargin(0.03);
  pad2.SetRightMargin(0.05);
  pad2.SetBottomMargin(0.20);

  pad1.cd();
  frame = obs->frame();
  auto frame2 = obs->frame();
  blind_data->plotOn(frame, RooFit::Name("data"));
  for(int ipdf = 0; ipdf < ncats; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    pdf->fitTo(*data, RooFit::PrintLevel(-1), RooFit::Warnings(0), RooFit::PrintEvalErrors(-1),
               RooFit::Range("LowSideband,HighSideband"), RooFit::NormRange("LowSideband,HighSideband"));
    pdf->plotOn(frame, RooFit::LineColor(colors[ipdf % ncolors]), RooFit::LineStyle((ipdf >= ncolors) ? kDashed : kSolid),
                RooFit::Title(titles[ipdf]), RooFit::Name(Form("pdf_%i", ipdf)), RooFit::Range("full"), RooFit::NormRange("LowSideband,HighSideband"));
    auto h = frame->pullHist("data", Form("pdf_%i", ipdf));
    h->SetLineWidth(2);
    h->SetLineColor(colors[ipdf % ncolors]);
    h->SetLineStyle((ipdf >= ncolors) ? kDashed : kSolid);
    h->SetMarkerSize(0.1);
    frame2->addPlotable(h, "PE1");
  }
  frame->Draw();
  leg->SetTextSize(0.036);
  leg->Draw("same");
  frame->GetYaxis()->SetRangeUser(0.1, 1.4*frame->GetMaximum());
  frame->GetXaxis()->SetLabelSize(0.);

  pad2.cd();
  frame2->Draw();

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame2->SetTitle("");
  frame2->GetYaxis()->SetTitle("#frac{N_{data} - N_{fit}}{#sigma_{data}}");

  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.2);

  frame2->GetYaxis()->SetRangeUser(-4., 4.);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->GetYaxis()->SetLabelSize(0.07);
  frame2->GetYaxis()->SetLabelOffset(0.01);
  frame2->GetYaxis()->SetTitleSize(0.07);
  frame2->GetYaxis()->SetTitleOffset(0.58);
  frame2->GetXaxis()->SetLabelSize(0.07);
  frame2->GetXaxis()->SetLabelOffset(0.008);
  frame2->GetXaxis()->SetTitleSize(0.07);
  frame2->GetXaxis()->SetTitleOffset(1.0);

  TLine line(frame->GetXaxis()->GetXmin(), 0., frame->GetXaxis()->GetXmax(), 0.);
  line.SetLineWidth(2);
  line.SetLineColor(kBlack);
  line.SetLineStyle(kDashed);
  line.Draw("same");
  pad2.SetGridy();

  c2.SaveAs(Form("inputs_refit_%s_%s.png", mp.Data(), bin));
  return 0;
}
