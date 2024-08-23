# Get the N(gen) information for the Z prime MC samples for normalization purposes
import os
import ROOT as rt

rt.gROOT.SetBatch(True)

### default path
path="/eos/cms/store/group/phys_smp/ZLFV/MC_generation/"

### MC signal mass points
sgn_masspoints=["200","400","600","800","1000"]
sgn_masspoint_files={
               "200" :"ZEMu_NANO_M200_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M200_2018_100k4/240710_144622/0000/*.root",\
               "400" :"ZEMu_NANO_M400_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M400_2018_100k4/240710_144531/0000/*.root",\
               "600" :"ZEMu_NANO_M600_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M600_2018_100k4/240710_144439/0000/*.root",\
               "800" :"ZEMu_NANO_M800_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M800_2018_100k4/240710_143659/0000/*.root",\
               "1000":"ZEMu_NANO_M1000_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M1000_2018_100k4/240710_143332/0000/*.root"}
norms = []
for mpoint in sgn_masspoints:
    # Create a TChain for the signal
    cc=rt.TChain("Events")
    filepath = path+sgn_masspoint_files[mpoint]
    print filepath
    cc.Add(filepath)
    norms.append(cc.GetEntries())
    print "Norm for mass point %-4s: %i" % (mpoint, norms[-1])

print norms
