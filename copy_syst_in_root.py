import ROOT as rt

nom_file="v4fit_shape_bin2-TH1.root"
nom_histos=["background","signal","data_obs"]
syst_fileUp="v4fit_shape_syst_Upbin2-TH1.root"
syst_fileDown="v4fit_shape_syst_Downbin2-TH1.root"
syst_histosUp=["background_fitUp"]
syst_histosDown=["background_fitDown"]
output_name="v4fit_shape_inclsys_bin2-TH1.root"


outhistos=[]

for hfile,histos in zip([nom_file,syst_fileUp,syst_fileDown],[nom_histos,syst_histosUp,syst_histosDown]):
  fin = rt.TFile(hfile,"READ")
  for histo in histos:
    hclone = (fin.Get(histo)).Clone()
    hclone.SetDirectory(0)
    hclone.SetName(histo)
    outhistos.append(hclone)
    #hclone.Write()
  fin.Close()

fout = rt.TFile(output_name,"RECREATE")
for histo in outhistos:
  histo.Write()
