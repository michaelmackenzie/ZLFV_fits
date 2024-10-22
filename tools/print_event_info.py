# Print event information in a TTree
import os
import argparse
import ROOT as rt

parser = argparse.ArgumentParser()
parser.add_argument("--input-file", dest="input_file",required=True, type=str, help="File to process")
parser.add_argument("--tree", dest="tree",default="mytreefit", type=str, help="Tree to process")
parser.add_argument("--max-entries", dest="max_entries",default=1, type=int, help="Maximum entries to print")
parser.add_argument("--first-entry", dest="first_entry",default=0, type=int, help="First entry process")
parser.add_argument("--cuts", dest="cuts",default="", type=str, help="Event cuts")

args, unknown = parser.parse_known_args()
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()

f = rt.TFile.Open(args.input_file, 'READ')
if not f: exit()
if not f.GetListOfKeys().Contains(args.tree):
    print 'Tree %s not found' % (args.tree)
    exit()
t = f.Get(args.tree)

t.Draw(">>elist", args.cuts, "goff")
elist = rt.gDirectory.Get('elist')

nentries = t.GetEntries()
nseen = 0
branches = [br.GetName() for br in t.GetListOfBranches()]

for entry in range(nentries):
    if not elist.Contains(entry): continue
    nseen += 1
    if nseen <= args.first_entry: continue
    print '>>> Entry %i (nseen = %i / %i)' % (entry, nseen-args.first_entry, args.max_entries)
    t.GetEntry(entry)
    for br in branches:
        print ' %-15s:' % (br)
        print getattr(t, br)
    if (nseen - args.first_entry) >= args.max_entries - 1: break
