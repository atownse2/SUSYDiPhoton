from ROOT import *
from utils import *
import os, sys, math
import numpy as np
from array import array
import pickle
gStyle.SetOptStat(0)
gROOT.SetBatch(1)

### Define data to be used
fDir = "hists/"

#dType = "DoubleEG"
#run = "2016"

#dType = "DoubleEG"
#run  = "2017"

#dType = "EGamma"
#run = "2018"

dType = "DYJetsToLL_M-50"
run = "Summer16v3"
rtag = "HT"


dName = dType + "_" + run

f_name = dName + "_fakerate.root"


### Define hist dicts and bins
from init_hists import *


##Load Saved Dicts
histDict = load_obj(dName + "_histDict")

for kVar in kVars:

  for tag in tags:

    bin_edges = varBins[kVar]
    if "eta" in kVar:
      var_bins = varBins[kVar]
      bin_edges = [0,1]
    else:
      var_bins = l2b(bin_edges)    

    h_title = dName + "_" + kVar + "_" + tag

    h_bins = array("d", bin_edges)
    h = TH1F(h_title, h_title, len(h_bins)-1, h_bins)

    if "eta" in kVar:
      h.GetXaxis().SetBinLabel(1, "barrel")
      #h.GetXaxis().SetBinLabel(2, "endcap")


    c = TCanvas(h_title+"_c", h_title+"_c", 800, 800)
    c.cd()
    
    h.GetXaxis().SetTitle("probe " + kVar)
    h.GetYaxis().SetTitle("fakerate")
    h.GetYaxis().SetRangeUser(0,0.06)
    h.SetMinimum(0)

    for i, b in enumerate(var_bins):
      h_ee_name = getHistName("ee", kVar, b, tag) 
      h_eg_name = getHistName("eg", kVar, b, tag)

      ee = histDict[h_ee_name]
      eg = histDict[h_eg_name]

      eeSig, eeSigError = ee["nSig"], ee["nSigError"]
      egSig, egSigError = eg["nSig"], ee["nSigError"]

      if eeSig*egSig == 0:
        print("Fit issue in " + kVar + " " + str(b))
        continue

      f, e = fakerate(eeSig, eeSigError, egSig, egSigError)

      if tag == "all" and kVar == "eta":
        print(b + " : fakerate: " + str(f) + " error: " + str(e))

      h.SetBinContent(i+1, f)
      h.SetBinError(i+1, e)

    h.Draw("histe")

    c.Update()

    c.SaveAs("plots/" + h_title + ".png")

for eVar in eVars:

  tag = None

  bin_edges = varBins[eVar]
  var_bins = l2b(bin_edges)    

  h_title = dName + "_" + eVar

  h_bins = array("d", bin_edges)
  h = TH1F(h_title, h_title, len(h_bins)-1, h_bins)

  c = TCanvas(h_title+"_c", h_title+"_c", 800, 800)
  c.cd()
    
  h.GetXaxis().SetTitle(h_title)
  h.GetYaxis().SetTitle("fakerate")
  h.GetYaxis().SetRangeUser(0,0.06)
  h.SetMinimum(0)

  for i, b in enumerate(var_bins):
    h_ee_name = getHistName("ee", eVar, b, tag) 
    h_eg_name = getHistName("eg", eVar, b, tag)

    ee = histDict[h_ee_name]
    eg = histDict[h_eg_name]

    eeSig, eeSigError = ee["nSig"], ee["nSigError"]
    egSig, egSigError = eg["nSig"], ee["nSigError"]

    if eeSig*egSig == 0:
      print("Fit issue in " + eVar + " " + str(b))
      continue

    f, e = fakerate(eeSig*2, eeSigError*2, egSig, egSigError) ##Weird factor of 2 needed to get same range of fakerates as the samples with kinematic variables, likely because the ee is filled for both electrons but met for example is only filled once per event

    h.SetBinContent(i+1, f)
    h.SetBinError(i+1, e)

  h.Draw()

  c.Update()

  c.SaveAs("plots/" + h_title + ".png")


"""
    line = TLine(sDict["bins"][0], fakerate, sDict["bins"][-1], fakerate )
    line.SetLineWidth(2)
    line.SetLineColor(kBlue)
    line.Draw()

    err = fakerate*0.67

    linep = TLine(sDict["bins"][0], fakerate+err, sDict["bins"][-1], fakerate+err )
    linep.SetLineWidth(2)
    linep.SetLineColor(kGreen)
    #linep.Draw()
    linem = TLine(sDict["bins"][0], fakerate-err, sDict["bins"][-1], fakerate-err )
    linem.SetLineWidth(2)
    linem.SetLineColor(kGreen)
    #linem.Draw()


    c.SaveAs("plots/" + title + ".png")
    c.Close()

    h.Write()
"""
