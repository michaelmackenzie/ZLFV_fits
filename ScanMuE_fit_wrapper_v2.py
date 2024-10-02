# Perform a scan of the emu data, fitting for a potential signal resonance
import os
import argparse
import ROOT as rt
from array import array
from ScanMuE_wrapper_helper import *
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')
from multiprocessing import Process

def eff_cor_wrt18(bdt_min,bdt_max,year):
    
    if bdt_min=="0.3" and bdt_max=="0.7" and year=="2016":
       return 1.018
    elif bdt_min=="0.3" and bdt_max=="0.7" and year=="2017":
       return 1.002
    elif bdt_min=="0.7" and bdt_max=="1.01" and year=="2016":
       return 0.9469 
    elif bdt_min=="0.7" and bdt_max=="1.01" and year=="2017":
       return 0.9837
    else: 
       return 1;



def proc_unit(sgn_argument, bkg_argument):
    if sgn_argument !='':
       os.system(sgn_argument)
    if bkg_argument != '':
       os.system(bkg_argument)


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="test", type=str,help="output root name")
parser.add_argument("--data-file", dest="data_file",default="", type=str,help="data file")
parser.add_argument("--xgb-min", dest="xgb_min",default="0.70", type=str,help="data file")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.01", type=str,help="data file")
parser.add_argument("--scan-min", dest="scan_min",default=110., type=float,help="mass signal")
parser.add_argument("--scan-max", dest="scan_max",default=500., type=float,help="mass signal")
parser.add_argument("--scan-step", dest="scan_step",default=1.0, type=float,help="mass signal")
parser.add_argument("--component", dest="component",type=str,default="all",help="Options: all sgn bkg ")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="data file")
parser.add_argument("--skip-shape-dc", dest="skip_shape_dc",default=False, action='store_true',help="shape experiment")
parser.add_argument("--year", dest="year",default="Run2", type=str,help="year")
parser.add_argument("--fit-version", dest="ver",default="2", type=str,help="fit version")
parser.add_argument("--outvar", dest="outvar",default="mass_ll",type=str,help="data file")
parser.add_argument("--param-name", dest="param_name",default="bin",type=str,help="data file")
parser.add_argument("--skip-fit", dest="skip_fit",default=False, action='store_true',help="fit skip")
parser.add_argument("--skip-sgn-syst", dest="skip_sgn_syst",default=False, action='store_true',help="shape experiment")
parser.add_argument("--skip-bkg-altfits", dest="skip_bkg_altfits",default=False, action='store_true',help="shape experiment")
parser.add_argument("--skip-dc", dest="skip_dc",default=False, action='store_true',help="Skip datacard creation")
parser.add_argument("--mass-point", dest="mass_point",default=-1,type=int,help="Single mass point to process")
parser.add_argument("--full-mass", dest="full_mass",default=False,action='store_true',help="Fit the entire mass distribution")
parser.add_argument("--log-files", dest="log_files",default=False,action='store_true',help="Write individual mass point fits to log files")
parser.add_argument("--dry-run", dest="dry_run",default=False,action='store_true',help="Don't execute script calls")

args, unknown = parser.parse_known_args()


############################# hard-set configuaration #########################

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exitting"
   exit()

MaxMasspoints=-1 #X for debug; -1 for run
sr_buffer = 1 #in signal width units
region_buffer = 10 #in signal width units

### MC signal mass points
sgn_masspoints=["100","125","150","175","200","300","400","500"]
sgn_masspoint_files={
               "100":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM100_mcRun18.root",\
               "125":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM125_mcRun18.root",\
               "150":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM150_mcRun18.root",\
               "175":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM175_mcRun18.root",\
               "200":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM200_mcRun18.root",\
               "300":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM300_mcRun18.root",\
               "400":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM400_mcRun18.root",\
               "500":"Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM500_mcRun18.root"}
ndens={"100":99200.,"125":98400.,"150":90500.,"175":89900.,"200":96300.,"300":99700.,"400":97600.,"500":98700.}



### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./FiguresScan/%s/" % (args.name)
carddir = "./DatacardsScan/%s/" % (args.name)

if os.path.exists(carddir):
   os.system("rm -rI "+carddir)

os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))
os.system("[ ! -d log ] && mkdir log")


#if not os.path.exists(carddir+"WorkspaceScanBKG"):
#   os.symlink("../../WorkspaceScanBKG", carddir+"WorkspaceScanBKG")
#if not os.path.exists(carddir+"WorkspaceScanSGN"):
#   os.symlink("../../WorkspaceScanSGN", carddir+"WorkspaceScanSGN")
   

### Data
data_files={"2016":"Meas_full_bdt_v7_emu_scan_data_Run16.root","2017":"Meas_full_bdt_v7_emu_scan_data_Run17.root","2018":"Meas_full_bdt_v7_emu_scan_data_Run18.root","Run2":"Meas_full_bdt_v7_emu_scan_data_Run1*.root"}
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":36.33*eff_cor_wrt18(args.xgb_min,args.xgb_max,"2016")+41.48*eff_cor_wrt18(args.xgb_min,args.xgb_max,"2017")+59.83}


### sf
#sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)"
sf="1"


### arange flags
unblind="false"
if args.unblind:
   unblind="true"
shape_dc="true"
if args.skip_shape_dc:
   shape_dc="false"
do_sgn_syst="true"
if args.skip_sgn_syst:
   do_sgn_syst="false"
do_bkg_altfits="true"
if args.skip_bkg_altfits:
   do_bkg_altfits="false"
if args.data_file=="":
   args.data_file=path+data_files[args.year]

rt.gROOT.SetBatch(True)

##############################################################################
################################### code #####################################

#----------------------------------------------
# Get the signal parameterization vs. mass
#----------------------------------------------
par_yield, par_width = mass_analysis(sgn_masspoints,path,sgn_masspoint_files, args.xgb_min, args.xgb_max, sf, ndens, lumis[args.year], figdir)


#----------------------------------------------
# Perform the signal scan at all mass points
#----------------------------------------------
###### main code ######
if MaxMasspoints>0:
   print "WILL PRODUCE",MaxMasspoints,"points"
NextPoint=True 

# set intitial scan mass
sr_center=args.scan_min

jobs = []

cnt=-1
script_head = 'echo ' if args.dry_run else ''

while (NextPoint):
  cnt+=1
  # calculate blind range and expected widths/yields for a mass
 # sr_width = round(par_width[0]+par_width[1]*sr_center+par_width[2]*sr_center*sr_center,2)
  sr_width = sr_center*1./50.
  sr_min = round(sr_center - sr_buffer*sr_width,2)
  sr_max = round(sr_center + sr_buffer*sr_width,2)
  sr_yld = round(par_yield[0] + par_yield[1]*sr_center + par_yield[2]*sr_center*sr_center,2)
  min_mass = round(sr_center - region_buffer*sr_width,2) if not args.full_mass else  95.
  max_mass = round(sr_center + region_buffer*sr_width,2) if not args.full_mass else 700.
  if min_mass<95: min_mass=95.
  
  print "SR #",cnt," central",sr_center,"width",sr_width,"min",sr_min,"max",sr_max,"yield",sr_yld

  if(args.mass_point < 0 or cnt == args.mass_point):        
    if not args.skip_fit:
       sgn_argument=''
       bkg_argument=''
       if args.component == "sgn" or args.component == "all":
          tail = (' >| log/fit_sgn_%s_mp%i.log' % (args.name, cnt)) if args.log_files else ''
          sgn_argument= script_head + 'root -l -b -q ScanMuE_fit_sgn_v'+args.ver+'.C\'("'+args.name+"_mp"+str(cnt)+'",'+str(min_mass)+','+str(max_mass)+','+str(sr_center)+','+str(sr_width)+','+str(sr_yld)+','+shape_dc+',"'+args.outvar+'",'+do_sgn_syst+',"'+args.param_name+'","'+carddir+'/WorkspaceSGN/")\'' + tail
       if args.component == "bkg" or args.component == "all":    
          tail = (' >| log/fit_bkg_%s_mp%i.log' % (args.name, cnt)) if args.log_files else ''
          bkg_argument= script_head + 'root -l -b -q ScanMuE_fit_bkg_v'+args.ver+'.C\'("'+args.name+"_mp"+str(cnt)+'","'+args.data_file+'","'+args.xgb_min+'","'+args.xgb_max+'",'+str(min_mass)+','+str(max_mass)+','+str(sr_min)+','+str(sr_max)+','+unblind+','+shape_dc+',"'+args.outvar+'","'+args.param_name+'","'+carddir+'/WorkspaceBKG/")\'' + tail
       job = Process(target = proc_unit, args=(sgn_argument,bkg_argument,))
       jobs.append(job)
                     
       # Create a corresponding datacard
    if not args.dry_run:
       print_datacard(carddir, args.name, args.ver, str(cnt), args.param_name, sr_center)

  # next iteration mass and exit conditions
  sr_center = round(sr_center +args.scan_step*sr_width,2)
  if cnt>= MaxMasspoints and MaxMasspoints>0:
     print "Requested only",MaxMasspoints,"run"
     NextPoint=False
  if sr_center >= args.scan_max:
     NextPoint=False

print("Parallel processing 8 threads for",len(jobs),"points")
for ithread in range(0,len(jobs),8):
    nthread = 8+ithread
    if (nthread>len(jobs)):
       nthread=len(jobs)-1
    print("Working", ithread,"to",nthread,'from',len(jobs),"done",int(ithread*100.0/len(jobs)),"%")   
    for job in jobs[ithread:nthread]:
      job.start()
    for job in jobs[ithread:nthread]:    
      job.join()
  
print "scannned",cnt,"points"



