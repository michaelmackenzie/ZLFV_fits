import ROOT as rt


Create = False
infiles =["template_zemu_embed_v3_bin1.root", "template_zemu_embed_v3_bin2.root", "template_zemu_embed_v3_bin3.root"]
outfiles =["bin1", "bin2", "bin3"]
Fzmm=0

Validate = True
embed= "template_zemu_embed_v3_bin3.root"
data= "/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root"
cuts_data = "0.9<xgb"
outname="bin3_logy"

if Create:
  for idx in range(len(infiles)):
  
    tfile = rt.TFile(infiles[idx],"READ")  
    stack = tfile.Get("bkg_stack")
    th1s = stack.GetHists()
    
    c1 = rt.TCanvas("c"+str(idx),"c1",800,600)
    leg = rt.TLegend(0.7,0.7,1,1)
    htotal = rt.TH1F()
    for ith1,th1 in enumerate(th1s):
      print th1.GetName()
      leg.AddEntry(th1,th1.GetName())
      th1.SetFillColor(0)
      th1.SetLineWidth(2)
      if ith1==0:
        htotal=th1.Clone()
        htotal.SetName("hbkg_safe")
        htotal.SetTitle("hbkg_safe")
        th1.Draw("HIST")
        th1.GetYaxis().SetRangeUser(0,1000)
      else:
        if "Z->" in th1.GetName():
          th1.Scale(Fzmm)
        htotal.Add(th1)
        th1.Draw("HIST sames")
    
    leg.Draw("sames")
    
    htotal.SetLineColor(1)
    htotal.Draw("HIST sames")
    leg.AddEntry(htotal,"total")
    
    leg.Draw("sames")
    c1.SaveAs("ctmpl_embed_ZmmF_"+str(Fzmm)+"_"+outfiles[idx]+".png")
    
    fout = rt.TFile("template_zemu_v3_custom_ZmmF_"+str(Fzmm)+"_"+outfiles[idx]+".root","RECREATE")
    htotal.Write();


if Validate:
   rt.gStyle.SetOptStat(0)
   tfile = rt.TFile(embed,"READ")
   hmc = tfile.Get("hbkg_safe")  
  
   cc_data = rt.TChain("mytreefit")
   cc_data.Add(data)
  
   hdata = rt.TH1F("hdata","",80,70,110) 
   cc_data.Draw("mass_ll>>hdata",cuts_data+" && ( ( 70<mass_ll && mass_ll<85) || ( 95<mass_ll &&  mass_ll<110) )");
   

   
   c1 = rt.TCanvas("c1","",800,600)
   pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
   pad1.SetBottomMargin(0)  # joins upper and lower plot
   pad1.SetGridx()
   pad1.Draw()
   # Lower ratio plot is pad2
   c1.cd()  # returns to main canvas before defining pad2
   pad2 = rt.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
   pad2.SetTopMargin(0)  # joins upper and lower plot
   pad2.SetBottomMargin(0.2)
   pad2.SetGridx()
   pad2.Draw()
   pad1.cd()
   pad1.SetLogy()
   hmc.Draw("HIST")
   hdata.SetLineWidth(2)
   hdata.Draw("PE sames")
   pad2.cd()
   hratio= hdata.Clone()
   hratio.SetName("hratio")
   hratio.Divide(hmc)
#   hratio.SetMarkerStyle(2)
   hratio.SetLineWidth(2)   
   hratio.GetYaxis().SetLabelSize(0.12)
   hratio.GetXaxis().SetLabelSize(0.15)
   hratio.GetYaxis().SetRangeUser(0.8,1.2)
   hratio.Draw("PE")
   c1.SaveAs("comp_data_embed_"+outname+".png")

