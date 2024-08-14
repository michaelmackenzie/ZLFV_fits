from ROOT import *


gROOT.SetBatch(True)
workfile = TFile("workspace_v15_bkg_bin1.root","READ")

ws = workfile.Get("workspace_background")

ws.Print()

datas = ws.allData()

lepm = ws.var("lepm_11")

xframe_data = lepm.frame()

leg_data = TLegend(0.7,0.7,1,1)
for idata,data in enumerate(datas):
  data.plotOn(xframe_data,RooFit.Name(data.GetName()), RooFit.MarkerColor(idata+1),RooFit.LineColor(idata+1))
  leg_data.AddEntry(xframe_data.findObject(data.GetName()),data.GetName(),"P")

c1 = TCanvas("c1","",800,600)
xframe_data.Draw("PE")
leg_data.Draw("sames")
c1.SaveAs("ws_view_datas.png")

xframe_pdf = lepm.frame()
data_pdf = ws.data("sim_binned_obs_bin1")
data_pdf.plotOn(xframe_data,RooFit.Name("data"),RooFit.Binning(80))

pdfs = ws.allPdfs()

leg_pdf = TLegend(0.7,0.7,1,1)
for ipdf,pdf in enumerate(pdfs):
  print pdf.GetName()
  if "bin1_" in pdf.GetName() or "bin2_" in pdf.GetName() or "bin3_" in pdf.GetName(): continue
  pdf.plotOn(xframe_pdf,RooFit.LineColor(ipdf+1),RooFit.Name(pdf.GetName()),RooFit.NormRange("70-110"),RooFit.Binning(80))
  leg_pdf.AddEntry(xframe_pdf.findObject(pdf.GetName()),pdf.GetName(),"L")


c2 = TCanvas("c2","",800,600)
xframe_pdf.Draw()
leg_pdf.Draw("sames")
c2.SaveAs("ws_view_pdfs.png")


vrs = ws.allVars()

for vr in vrs:
  vr.Print()


