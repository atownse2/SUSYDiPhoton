from ROOT import * 
import sys
import numpy as np
from array import array
from os import listdir
from os.path import isfile, join

from collections import OrderedDict

gStyle.SetOptStat(0)
gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

fdir = "hists/"


#dType = "DYJetsToLL_M-50"
#run = "Summer16v3"

dType = "WGJets"
run = "Summer16v3"

fin_name = fdir + dType + "_" + run + "_genfakerate_ns.root"
fout_name = fdir + dType + "_" + run + "_genfakerate.root"


fin = TFile.Open(fin_name) 

if not fin:
  print("File not opened")
  quit()


### Define hist dicts and bins
ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets_DiLept_TuneCUETP8M1": ["_"]}

varBins = { "pt"    : [80,100,120,150,200,300],
            "eta"   : np.linspace(-2.4, 2.4, num=24),
            "njets"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250]}

genTypes = ["genEle"]

recoTypes = ["recoPho", "recoEle"]

regions = ["barrel" , "endcap"]

eVars = ["pt", "eta" , "njets" , "met"]

ptRanges = ptRangeDict[dType]

def getHistName(genType, recoType, region, var):
  return genType + "_" + recoType + "_" + region + "_" + var

hists = OrderedDict()

for genType in genTypes:
  for recoType in recoTypes:
    for region in regions:
      for eVar in eVars:
        h_name = getHistName(genType, recoType, region, eVar)
        h_bins = array("d", varBins[eVar])
        hists[h_name] = TH1F( h_name, h_name, len(h_bins)-1, h_bins)

hEventsDict = {ptRange: fin.Get("hEvents_"+ptRange) for ptRange in ptRanges}

for h_name, hist in hists.items():
  for ptRange, hEvents, in hEventsDict.items():
      hName_pt = h_name + "_" + ptRange
      h = fin.Get(hName_pt)
      h.Scale(2.5/hEvents.GetEntries())
      hist.Add(h)

### Write Hists
fout = TFile.Open(fout_name, "RECREATE")
fout.cd()

for hName, hist in hists.items():
  hist.Write()

### Make plots
for genType in genTypes:
  for region in regions:
    for eVar in eVars:
      title = dType + "_" + genType + "_" + region + "_" + eVar

      h_recoPho = hists[getHistName(genType, "recoPho" , region, eVar)]
      h_recoEle = hists[getHistName(genType, "recoEle" , region, eVar)]

      c = TCanvas(title+"_c", title+"_c", 1000, 1000)
      c.SetLogy()
      l = TLegend(0.6,0.8,0.9,0.9)
      c.cd()

      h_recoPho.SetStats(0)
      h_recoEle.SetStats(0)

      h_recoEle.SetLineColor(2)

      l.AddEntry(h_recoPho, "recoPho")
      l.AddEntry(h_recoEle, "recoEle")

      recoPhoMax = h_recoPho.GetMaximum() + h_recoPho.GetBinError(h_recoPho.GetMaximumBin())
      recoEleMax = h_recoEle.GetMaximum() + h_recoEle.GetBinError(h_recoEle.GetMaximumBin())

      yMax = max(recoPhoMax,recoEleMax)*1.0

      hRatio = TRatioPlot(h_recoPho,h_recoEle) 
      hRatio.Draw()

      hRatio.GetUpperPad().cd()

      h_recoPho.Draw("hist E1 same")
      h_recoEle.Draw("hist E1 same")

      hRatio.GetUpperRefYaxis().SetRangeUser(1,yMax)
      hRatio.GetLowerRefGraph().SetMaximum(0.025)

      

      c.Update()
      c.SaveAs("plots/" + title + ".png")


fin.Close()
fout.Close()
