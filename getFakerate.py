from ROOT import *
from utils import *
import os, sys, math
import numpy as np

gStyle.SetOptStat(0)
gROOT.SetBatch(1)

fDir = "hists/"

dType = "DYJetsToLL_M-50"
run = "Summer16v3"

#dType = "DoubleEG"
#run = "2016"

#dType = "DoubleEG"
#run  = "2017"

#dType = "EGamma"
#run = "2018"

#fName = fDir + dType + "_Run" + run + "_fakerate.root" 

fName = dType + "_" + run +  "_fakerate.root"

f = TFile(fDir + fName)

selections = ["ee", "eg"]

tags = ["tagHpt", "tagLpt"]

regions = ["barrel", "endcap"]

dataDict = { }

for region in regions:
  for selection in selections:
    hName = selection + "invmass_" + region
    dataDict[hName] = {}
    h = TH1F(hName, hName, 128, 0, 150)
    for tag in tags:
      hName_tag = selection + "invmass_" + tag + "_eta_" + region
      h_tmp = f.Get(hName_tag)
      h_tmp.Sumw2()
      h.Add(h_tmp,1) 
    dataDict[hName]["hist"] = h


for hName, hDict in dataDict.items():
  cName = dType + "_Run" + run + "_" + hName
  c = TCanvas("c", cName, 600, 600)
  
  c.cd()

  nSig, nSigError, fitStatus, frame = fitHist(hDict["hist"])

  hDict["nSig"] = nSig
  hDict["nSigError"] = nSigError

  frame.Draw()

  label = TPaveText( 0.6, 0.3, 0.9, 0.80, 'NDC') 
  label.AddText("N Signal = " + str(nSig))
  label.AddText("Err Sig  = " + str(nSigError))
  label.AddText("Fit Status: " + str(fitStatus))
  label.Draw()

  c.SaveAs("fits/" + cName + ".png")
  c.Close()


##Calculate the fakerate
for region in regions:

  eeDict = dataDict["eeinvmass_" + region]
  egDict = dataDict["eginvmass_" + region]

  eeSig = eeDict["nSig"]
  eeSigError = eeDict["nSigError"]

  egSig = egDict["nSig"]
  egSigError = egDict["nSigError"]

  f, e = fakerate(eeSig, eeSigError, egSig, egSigError)

  print(region + " fakerate: " + str(f) + " error: " + str(e))
