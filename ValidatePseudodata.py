import ROOT as rt
from tdrstyle import *
from cms_lumi import *


setTDRStyle()

def create_canvas(name):
  H_ref = 600;
  W_ref = 800;
  W = W_ref
  H  = H_ref
  T = 0.08*H_ref
  B = 0.12*H_ref
  L = 0.12*W_ref
  R = 0.04*W_ref
  canvas = rt.TCanvas(name,name,50,50,W,H)
  canvas.SetFillColor(0)
  canvas.SetBorderMode(0)
  canvas.SetFrameFillStyle(0)
  canvas.SetFrameBorderMode(0)
  canvas.SetLeftMargin( L/W )
  canvas.SetRightMargin( R/W )
  canvas.SetTopMargin( T/H )
  canvas.SetBottomMargin( B/H )
  canvas.SetTickx(0)
  canvas.SetTicky(0)
  return canvas


def load_histo(histo,color):
    histo.GetXaxis().SetNdivisions(6,5,0)
    histo.GetYaxis().SetNdivisions(6,5,0)
    histo.GetXaxis().SetTitleOffset(0.95)
    histo.GetYaxis().SetTitleOffset(0.95)
    histo.SetLineColor(color)
    histo.SetMarkerColor(color)
    histo.SetLineWidth(3)

def plot_histos(histos,labels,colors,name,xaxis,norm=False, logy=False, legPos="TR",rng=[],yaxis=None):
   rt.gStyle.SetOptStat(0)
   c1 =create_canvas(name)
   if legPos=="TR":
      leg = rt.TLegend(0.6,0.7,0.9,0.9)
   elif legPos=="TL":
      leg = rt.TLegend(0.2,0.7,0.4,0.9)
   for ihst in range(len(histos)):
     histo = histos[ihst]
     histo.SetMinimum(0)
     load_histo(histo,colors[ihst])
     if norm:
        histo.Scale(1./histo.Integral())

     if len(rng)>0:
        histo.GetYaxis().SetRangeUser(rng[0],rng[1])

     if ihst==0:
        histo.Draw("HIST")
        histo.GetXaxis().SetTitle(xaxis)
        if not yaxis == None:
           histo.GetYaxis().SetTitle("Events / "+yaxis+" GeV")
        else:
           histo.GetYaxis().SetTitle("Events")
     else:
        histo.Draw("HIST sames")

     if len(histos)>1:
        leg.AddEntry(histo,labels[ihst])

   if len(histos)>1:
     leg.Draw("sames")
   if logy:
     c1.SetLogy()

   CMS_lumi(c1, 5,  0 , aLittleExtra=0.07)
   c1.SaveAs(name+".png")



rt.gROOT.SetBatch(True)
tree_pseudo=rt.TChain("mytreefit")
tree_pseudo.Add("pseudo_data_from_MC_v2_r0.root")
tree_data=rt.TChain("mytreefit")
tree_data.Add("/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root")

hmll_psd = rt.TH1F("mll_psd","",40,70,110)
hmll_dt = rt.TH1F("mll_dt","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_psd","NormGen_wt*(0.3<xgb && xgb<1.0 && (mass_ll<85 || mass_ll>95))")
tree_data.Draw("mass_ll>>mll_dt","0.3<xgb && xgb<1.0 && (mass_ll<85 || mass_ll>95)")
plot_histos([hmll_psd,hmll_dt],["MC","Data"],[1,2],"psd_val_incl","M(e,#mu)",True)



hmll_bin1_psd = rt.TH1F("mll_bin1_psd","",40,70,110)
hmll_bin1_dt = rt.TH1F("mll_bin1_dt","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_bin1_psd","NormGen_wt*(0.3<xgb && xgb<0.7 && (mass_ll<85 || mass_ll>95))")
tree_data.Draw("mass_ll>>mll_bin1_dt","0.3<xgb && xgb<0.7 && (mass_ll<85 || mass_ll>95)")
plot_histos([hmll_bin1_psd,hmll_bin1_dt],["MC","Data"],[1,2],"psd_val_bin1","M(e,#mu)",True)


hmll_bin2_psd = rt.TH1F("mll_bin2_psd","",40,70,110)
hmll_bin2_dt = rt.TH1F("mll_bin2_dt","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_bin2_psd","NormGen_wt*(0.7<xgb && xgb<0.9 && (mass_ll<85 || mass_ll>95))")
tree_data.Draw("mass_ll>>mll_bin2_dt","0.7<xgb && xgb<0.9 && (mass_ll<85 || mass_ll>95)")
plot_histos([hmll_bin2_psd,hmll_bin2_dt],["MC","Data"],[1,2],"psd_val_bin2","M(e,#mu)",True)


hmll_bin3_psd = rt.TH1F("mll_bin3_psd","",40,70,110)
hmll_bin3_dt = rt.TH1F("mll_bin3_dt","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_bin3_psd","NormGen_wt*(0.9<xgb && xgb<1.01 && (mass_ll<85 || mass_ll>95))")
tree_data.Draw("mass_ll>>mll_bin3_dt","0.9<xgb && xgb<1.01 && (mass_ll<85 || mass_ll>95)")
plot_histos([hmll_bin3_psd,hmll_bin3_dt],["MC","Data"],[1,2],"psd_val_bin3","M(e,#mu)",True)


hmll_tot = rt.TH1F("mll_tot","",40,70,110)
hmll_zmm = rt.TH1F("mll_zmm","",40,70,110)
hmll_zee = rt.TH1F("mll_zee","",40,70,110)
hmll_ztt = rt.TH1F("mll_ztt","",40,70,110)
hmll_nody = rt.TH1F("mll_nody","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_tot","NormGen_wt*(0.3<xgb && xgb<1.01)")
tree_pseudo.Draw("mass_ll>>mll_zmm","NormGen_wt*(0.3<xgb && xgb<1.01 && IsZmm==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_zee","NormGen_wt*(0.3<xgb && xgb<1.01 && IsZee==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_ztt","NormGen_wt*(0.3<xgb && xgb<1.01 && IsZtt==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_nody","NormGen_wt*(0.3<xgb && xgb<1.01 && dy==0)")
plot_histos([hmll_tot,hmll_zmm,hmll_zee,hmll_ztt,hmll_nody],["Total","Z->#mu#mu","Z->ee","Z->#tau#tau","WW/ttbar"],[1,2,3,4,5],"fakes_breakdown_incl","M(e,#mu)",False)

hmll_bin1_tot = rt.TH1F("mll_bin1_tot","",40,70,110)
hmll_bin1_zmm = rt.TH1F("mll_bin1_zmm","",40,70,110)
hmll_bin1_zee = rt.TH1F("mll_bin1_zee","",40,70,110)
hmll_bin1_ztt = rt.TH1F("mll_bin1_ztt","",40,70,110)
hmll_bin1_nody = rt.TH1F("mll_bin1_nody","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_bin1_tot","NormGen_wt*(0.3<xgb && xgb<0.7)")
tree_pseudo.Draw("mass_ll>>mll_bin1_zmm","NormGen_wt*(0.3<xgb && xgb<0.7 && IsZmm==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin1_zee","NormGen_wt*(0.3<xgb && xgb<0.7 && IsZee==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin1_ztt","NormGen_wt*(0.3<xgb && xgb<0.7 && IsZtt==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin1_nody","NormGen_wt*(0.3<xgb && xgb<0.7 && dy==0)")
plot_histos([hmll_bin1_tot,hmll_bin1_zmm,hmll_bin1_zee,hmll_bin1_ztt,hmll_bin1_nody],["Total","Z->#mu#mu","Z->ee","Z->#tau#tau","WW/ttbar"],[1,2,3,4,5],"fakes_breakdown_bin1","M(e,#mu)",False)


hmll_bin2_tot = rt.TH1F("mll_bin2_tot","",40,70,110)
hmll_bin2_zmm = rt.TH1F("mll_bin2_zmm","",40,70,110)
hmll_bin2_zee = rt.TH1F("mll_bin2_zee","",40,70,110)
hmll_bin2_ztt = rt.TH1F("mll_bin2_ztt","",40,70,110)
hmll_bin2_nody = rt.TH1F("mll_bin2_nody","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_bin2_tot","NormGen_wt*(0.7<xgb && xgb<0.9)")
tree_pseudo.Draw("mass_ll>>mll_bin2_zmm","NormGen_wt*(0.7<xgb && xgb<0.9 && IsZmm==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin2_zee","NormGen_wt*(0.7<xgb && xgb<0.9 && IsZee==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin2_ztt","NormGen_wt*(0.7<xgb && xgb<0.9 && IsZtt==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin2_nody","NormGen_wt*(0.7<xgb && xgb<0.9 && dy==0)")
plot_histos([hmll_bin2_tot,hmll_bin2_zmm,hmll_bin2_zee,hmll_bin2_ztt,hmll_bin2_nody],["Total","Z->#mu#mu","Z->ee","Z->#tau#tau","WW/ttbar"],[1,2,3,4,5],"fakes_breakdown_bin2","M(e,#mu)",False)


hmll_bin3_tot = rt.TH1F("mll_bin3_tot","",40,70,110)
hmll_bin3_zmm = rt.TH1F("mll_bin3_zmm","",40,70,110)
hmll_bin3_zee = rt.TH1F("mll_bin3_zee","",40,70,110)
hmll_bin3_ztt = rt.TH1F("mll_bin3_ztt","",40,70,110)
hmll_bin3_nody = rt.TH1F("mll_bin3_nody","",40,70,110)
tree_pseudo.Draw("mass_ll>>mll_bin3_tot","NormGen_wt*(0.9<xgb && xgb<1.01)")
tree_pseudo.Draw("mass_ll>>mll_bin3_zmm","NormGen_wt*(0.9<xgb && xgb<1.01 && IsZmm==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin3_zee","NormGen_wt*(0.9<xgb && xgb<1.01 && IsZee==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin3_ztt","NormGen_wt*(0.9<xgb && xgb<1.01 && IsZtt==1 && dy==1)")
tree_pseudo.Draw("mass_ll>>mll_bin3_nody","NormGen_wt*(0.9<xgb && xgb<1.01 && dy==0)")
plot_histos([hmll_bin3_tot,hmll_bin3_zmm,hmll_bin3_zee,hmll_bin3_ztt,hmll_bin3_nody],["Total","Z->#mu#mu","Z->ee","Z->#tau#tau","WW/ttbar"],[1,2,3,4,5],"fakes_breakdown_bin3","M(e,#mu)",False)



tree_pseudo_fakes2 = rt.TChain("mytreefit")
tree_pseudo_fakes2.Add("pseudo_data_from_MC_v2_r0_ZmmR1.5.root")
tree_pseudo_fakes0p5 = rt.TChain("mytreefit")
tree_pseudo_fakes0p5.Add("pseudo_data_from_MC_v2_r0_ZmmR0.5.root")
hmll_bin1_fakes2_zmm = rt.TH1F("mll_bin1_fakes2_zmm","",40,70,110)
hmll_bin2_fakes2_zmm = rt.TH1F("mll_bin2_fakes2_zmm","",40,70,110)
hmll_bin3_fakes2_zmm = rt.TH1F("mll_bin3_fakes2_zmm","",40,70,110)
hmll_bin1_fakes2_tot = rt.TH1F("mll_bin1_fakes2_tot","",40,70,110)
hmll_bin2_fakes2_tot = rt.TH1F("mll_bin2_fakes2_tot","",40,70,110)
hmll_bin3_fakes2_tot = rt.TH1F("mll_bin3_fakes2_tot","",40,70,110)

hmll_bin1_fakes0p5_zmm = rt.TH1F("mll_bin1_fakes0p5_zmm","",40,70,110)
hmll_bin2_fakes0p5_zmm = rt.TH1F("mll_bin2_fakes0p5_zmm","",40,70,110)
hmll_bin3_fakes0p5_zmm = rt.TH1F("mll_bin3_fakes0p5_zmm","",40,70,110)
hmll_bin1_fakes0p5_tot = rt.TH1F("mll_bin1_fakes0p5_tot","",40,70,110)
hmll_bin2_fakes0p5_tot = rt.TH1F("mll_bin2_fakes0p5_tot","",40,70,110)
hmll_bin3_fakes0p5_tot = rt.TH1F("mll_bin3_fakes0p5_tot","",40,70,110)

tree_pseudo_fakes2.Draw("mass_ll>>mll_bin1_fakes2_zmm","NormGen_wt*(0.3<xgb && xgb<0.7 && IsZmm==1 && dy==1)")
tree_pseudo_fakes2.Draw("mass_ll>>mll_bin2_fakes2_zmm","NormGen_wt*(0.7<xgb && xgb<0.9 && IsZmm==1 && dy==1)")
tree_pseudo_fakes2.Draw("mass_ll>>mll_bin3_fakes2_zmm","NormGen_wt*(0.9<xgb && xgb<1.01 && IsZmm==1 && dy==1)")
tree_pseudo_fakes2.Draw("mass_ll>>mll_bin1_fakes2_tot","NormGen_wt*(0.3<xgb && xgb<0.7)")
tree_pseudo_fakes2.Draw("mass_ll>>mll_bin2_fakes2_tot","NormGen_wt*(0.7<xgb && xgb<0.9)")
tree_pseudo_fakes2.Draw("mass_ll>>mll_bin3_fakes2_tot","NormGen_wt*(0.9<xgb && xgb<1.01)")
tree_pseudo_fakes0p5.Draw("mass_ll>>mll_bin1_fakes0p5_zmm","NormGen_wt*(0.3<xgb && xgb<0.7 && IsZmm==1 && dy==1)")
tree_pseudo_fakes0p5.Draw("mass_ll>>mll_bin2_fakes0p5_zmm","NormGen_wt*(0.7<xgb && xgb<0.9 && IsZmm==1 && dy==1)")
tree_pseudo_fakes0p5.Draw("mass_ll>>mll_bin3_fakes0p5_zmm","NormGen_wt*(0.9<xgb && xgb<1.01 && IsZmm==1 && dy==1)")
tree_pseudo_fakes0p5.Draw("mass_ll>>mll_bin1_fakes0p5_tot","NormGen_wt*(0.3<xgb && xgb<0.7)")
tree_pseudo_fakes0p5.Draw("mass_ll>>mll_bin2_fakes0p5_tot","NormGen_wt*(0.7<xgb && xgb<0.9)")
tree_pseudo_fakes0p5.Draw("mass_ll>>mll_bin3_fakes0p5_tot","NormGen_wt*(0.9<xgb && xgb<1.01)")


plot_histos([hmll_bin1_tot,hmll_bin1_zmm,hmll_bin1_fakes2_tot,hmll_bin1_fakes2_zmm],["Total(Z->#mu#mu x1)","Z->#mu#mu(x1)","Total(Z->#mu#mu x1.5)","Z->#mu#mu(x1.5)"],[1,2,3,4],"double_fakes_bin1","M(e,#mu)",False)
plot_histos([hmll_bin2_tot,hmll_bin2_zmm,hmll_bin2_fakes2_tot,hmll_bin2_fakes2_zmm],["Total(Z->#mu#mu x1)","Z->#mu#mu(x1)","Total(Z->#mu#mu x1.5)","Z->#mu#mu(x1.5)"],[1,2,3,4],"double_fakes_bin2","M(e,#mu)",False)
plot_histos([hmll_bin3_tot,hmll_bin3_zmm,hmll_bin3_fakes2_tot,hmll_bin3_fakes2_zmm],["Total(Z->#mu#mu x1)","Z->#mu#mu(x1)","Total(Z->#mu#mu x1.5)","Z->#mu#mu(x1.5)"],[1,2,3,4],"double_fakes_bin3","M(e,#mu)",False)


plot_histos([hmll_bin1_tot,hmll_bin1_zmm,hmll_bin1_fakes0p5_tot,hmll_bin1_fakes0p5_zmm],["Total(Z->#mu#mu x1)","Z->#mu#mu(x1)","Total(Z->#mu#mu /2)","Z->#mu#mu(/2)"],[1,2,3,4],"half_fakes_bin1","M(e,#mu)",False)
plot_histos([hmll_bin2_tot,hmll_bin2_zmm,hmll_bin2_fakes0p5_tot,hmll_bin2_fakes0p5_zmm],["Total(Z->#mu#mu x1)","Z->#mu#mu(x1)","Total(Z->#mu#mu /2)","Z->#mu#mu(/2)"],[1,2,3,4],"half_fakes_bin2","M(e,#mu)",False)
plot_histos([hmll_bin3_tot,hmll_bin3_zmm,hmll_bin3_fakes0p5_tot,hmll_bin3_fakes0p5_zmm],["Total(Z->#mu#mu x1)","Z->#mu#mu(x1)","Total(Z->#mu#mu /2)","Z->#mu#mu(/2)"],[1,2,3,4],"half_fakes_bin3","M(e,#mu)",False)


tree_newid=rt.TChain("mytreefit")
tree_newid.Add("pseudo_data_from_MC_v2_r0_updateID.root")

hmll_bin1_newid_zmm = rt.TH1F("mll_bin1_newid_zmm","",40,70,110)
hmll_bin2_newid_zmm = rt.TH1F("mll_bin2_newid_zmm","",40,70,110)
hmll_bin3_newid_zmm = rt.TH1F("mll_bin3_newid_zmm","",40,70,110)
hmll_bin1_newid_tot = rt.TH1F("mll_bin1_newid_tot","",40,70,110)
hmll_bin2_newid_tot = rt.TH1F("mll_bin2_newid_tot","",40,70,110)
hmll_bin3_newid_tot = rt.TH1F("mll_bin3_newid_tot","",40,70,110)

tree_newid.Draw("mass_ll>>mll_bin1_newid_zmm","NormGen_wt*(0.3<xgb && xgb<0.7 && IsZmm==1 && dy==1)")
tree_newid.Draw("mass_ll>>mll_bin2_newid_zmm","NormGen_wt*(0.7<xgb && xgb<0.9 && IsZmm==1 && dy==1)")
tree_newid.Draw("mass_ll>>mll_bin3_newid_zmm","NormGen_wt*(0.9<xgb && xgb<1.01 && IsZmm==1 && dy==1)")
tree_newid.Draw("mass_ll>>mll_bin1_newid_tot","NormGen_wt*(0.3<xgb && xgb<0.7)")
tree_newid.Draw("mass_ll>>mll_bin2_newid_tot","NormGen_wt*(0.7<xgb && xgb<0.9)")
tree_newid.Draw("mass_ll>>mll_bin3_newid_tot","NormGen_wt*(0.9<xgb && xgb<1.01)")

plot_histos([hmll_bin1_tot,hmll_bin1_zmm,hmll_bin1_newid_tot,hmll_bin1_newid_zmm],["Total(old ID)","Z->#mu#mu(old ID)","Total(new ID)","Z->#mu#mu(new ID)"],[1,2,3,4],"half_newid_bin1","M(e,#mu)",False)
plot_histos([hmll_bin2_tot,hmll_bin2_zmm,hmll_bin2_newid_tot,hmll_bin2_newid_zmm],["Total(old ID)","Z->#mu#mu(old ID)","Total(new ID)","Z->#mu#mu(new ID)"],[1,2,3,4],"half_newid_bin2","M(e,#mu)",False)
plot_histos([hmll_bin3_tot,hmll_bin3_zmm,hmll_bin3_newid_tot,hmll_bin3_newid_zmm],["Total(old ID)","Z->#mu#mu(old ID)","Total(new ID)","Z->#mu#mu(new ID)"],[1,2,3,4],"half_newid_bin3","M(e,#mu)",False)
print "bin1: Z->mm reduction",1-hmll_bin1_newid_zmm.Integral()/hmll_bin1_zmm.Integral(),"total bkg",1-hmll_bin1_newid_tot.Integral()/hmll_bin1_tot.Integral()
print "bin2: Z->mm reduction",1-hmll_bin2_newid_zmm.Integral()/hmll_bin2_zmm.Integral(),"total bkg",1-hmll_bin2_newid_tot.Integral()/hmll_bin2_tot.Integral()
print "bin3: Z->mm reduction",1-hmll_bin3_newid_zmm.Integral()/hmll_bin3_zmm.Integral(),"total bkg",1-hmll_bin3_newid_tot.Integral()/hmll_bin3_tot.Integral()

