from ROOT import gStyle, gROOT, TCanvas, gPad
from utils import get_file_name, load_hists, genTypes, recoTypes, regions2, finalStates, eVars

import pandas as pd

gStyle.SetOptStat(0)
gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

fdir = "hists/"


#dType = "DYJetsToLL_M-50"
#era = "Summer16v3"

dType = "WGJets"
era = "Summer16v3"

#dType = "TTJets"
#era = "Summer16v3"

#ntuple_version = "TreeMaker"
ntuple_version = "TreeMakerRandS_skimsv8"

version = "060922v1"
#version_tag = None

fin_name = fdir + get_file_name(dType, era, ntuple_version, "fakeBKG", version)

print("File In: " + fin_name)


# Load hists
hists = load_hists(dType, fin_name)

for region in regions2:
  for fstate in finalStates:
    h_name = fstate + "_" + region + "_matrix"
    h = hists[h_name]

    c = TCanvas(h_name + "_c", h_name + "_c", 1000, 800)
    c.cd()

    x_title = "lead"
    y_title = "trail"

    x_title += "_recoEle" if fstate[0] == "e" else "_recoPho"
    y_title += "_recoEle" if fstate[1] == "e" else "_recoPho"

    x_title += "_barrel" if region[0] == "b" else "_endcap"
    y_title += "_barrel" if region[1] == "b" else "_endcap"

    h.GetXaxis().LabelsOption("a")
    h.GetYaxis().LabelsOption("a")

    h.GetXaxis().SetTitle(x_title)
    h.GetYaxis().SetTitle(y_title)

    h.GetYaxis().SetTitleOffset(2.0)
    gPad.SetLeftMargin(0.15)

    h.SetStats(False)
    h.Draw("TEXT")

    fName = dType + "_" + era + "_" + ntuple_version + "_" + version + "_" + h_name
    h.SetTitle(fName)
 
    c.Update()
    c.SaveAs("plots/" +fName+".png")


"""
for genType in genTypes:
  for region in regions:
    h_nEle = hists[getHistName(genType, "recoEle" , region, "nEle")]
    h_nPho = hists[getHistName(genType, "recoPho" , region, "nPho")]
    tf = h_nPho.GetEntries()/h_nEle.GetEntries()

    tfError = ratioPoissonError(h_nPho.GetEntries(), h_nEle.GetEntries())


    print("Transfer facter = nPho/nEle = " + str((tf, tfError)) + " " + region)


### Make plots
for genType in genTypes:
  for region in regions:

    for eVar in eVars:
      if eVar == "nEle" or eVar == "nPho":
        continue

      title = dType + "_" + era + "_" + ntuple_version + "_"  + genType + "_" + region + "_" + eVar
  
      if version_tag is not None:
        title +=  "_" + version_tag
 
      h_recoPho = hists[getHistName(genType, "recoPho" , region, eVar)]
      h_recoEle = hists[getHistName(genType, "recoEle" , region, eVar)]

      c = TCanvas(title+"_c", title+"_c", 1000, 800)
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

      yMax = max(recoPhoMax,recoEleMax)*1.25

      hRatio = TRatioPlot(h_recoPho,h_recoEle) 
      hRatio.Draw()

      l.Draw()

      hRatio.GetUpperPad().cd()

      h_recoPho.Draw("hist E1 same")
      h_recoEle.Draw("hist E1 same")

      hRatio.GetUpperRefYaxis().SetRangeUser(1,yMax)
      hRatio.GetLowerRefGraph().SetMinimum(0)
      hRatio.GetLowerRefGraph().SetMaximum(0.05)

 

      c.Update()
      c.SaveAs("plots/" + title + ".png")
"""

