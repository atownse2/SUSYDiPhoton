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
    for tag in tags:
      hName = selection + "invmass_" + tag + "_eta_" + region
      dataDict[hName] = {} 
      h = f.Get(hName)
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

  for tag in tags:

    eeDict = dataDict["eeinvmass_" + tag + "_eta_" + region]
    egDict = dataDict["eginvmass_" + tag + "_eta_" + region]

    eeSig = eeDict["nSig"]
    eeSigError = eeDict["nSigError"]

    egSig = egDict["nSig"]
    egSigError = egDict["nSigError"]

    f, e = fakerate(eeSig, eeSigError, egSig, egSigError)

    print(region + " " + tag + " fakerate: " + str(f) + " error: " + str(e))


