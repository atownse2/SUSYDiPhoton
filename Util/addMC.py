from ROOT import TCanvas, TFile, TH1F, PyConfig, TLegend, TRatioPlot
from ROOT import gROOT
import sys
import numpy as np
from array import array
from os import listdir
from os.path import isfile, join

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

fdir = "hists/"
dType = "DYJetsToLL_M-50"
run = "Summer16v3"

fin_name = fdir + dType + "_" + run + "_fakerate_ns.root"
fout_name = fdir + dType + "_" + run + "_fakerate.root"


fin = TFile.Open(fin_name) 

if not fin:
  print("File not opened")
  quit()

###

### Initialize hists
from init_hists import *
import init_hists_mc as mc

init_hists()

ptRanges = mc.ptRangeDict[dType]

hEventsDict = {ptRange: fin.Get("hEvents_"+ptRange) for ptRange in ptRanges}

for hName, hDict in histDict.items():
  for ptRange, hEvents, in hEventsDict.items():
      hName_pt = hName.replace("invmass_", "invmass_"+ptRange+"_")

      h = fin.Get(hName_pt)

      h.Scale(2.5/hEvents.GetEntries())

      hDict["hist"].Add(h)

#fin.Close() ##This breaks the code for some reason

### Write Hists
fout = TFile.Open(fout_name, "RECREATE")
fout.cd()

for hName, hDict in histDict.items():
  hDict["hist"].Write()

fin.Close()
fout.Close()
