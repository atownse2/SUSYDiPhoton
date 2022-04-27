from ROOT import *
from utils import *
import os, sys, math
import numpy as np

gStyle.SetOptStat(0)
gROOT.SetBatch(1)

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

fName = dName + "_fakerate.root"

#Create hist dict structure and load hists
from init_hists import *
#load_hists(fDir+fName)



infile = TFile(fDir+fName)

### Loop over hists and fit/plot
t_hists = ["eeinvmass", "eginvmass", "gginvmass"]



for hName in t_hists:
  cName = dName + "_" + hName
  c = TCanvas(cName, cName, 400, 400)

  c.cd()

  hist = infile.Get(hName) 

  nSig, nSigError, fitStatus, frame = fitHist(hist)

  histDict[hName]["nSig"] = nSig
  histDict[hName]["nSigError"] = nSigError

  frame.Draw()

  #label = TPaveText( 0.9, 0.7, 0.9, 0.80, 'NDC')
  label = TPaveText( 0.5, 0.5, 0.9, 0.9, "NDC" )  
  label.AddText("N Signal = " + str(nSig))
  label.AddText("Err Sig  = " + str(nSigError))
  label.AddText("Fit Status: " + str(fitStatus))
  label.Draw()

  c.SaveAs("fits/" + cName + ".png")
  c.Close()



for kVar in kVars:
  for tag in tags:
    for eSel in eventSelections:
      h_tag = eSel + "invmass" + "_" + tag + "_" + kVar
      hNames = [h for h in histDict.keys() if h_tag in h]

      cName = dName + "_" + h_tag
      c = TCanvas(cName, cName, 900, 900)
      
      c.cd()
      c.Divide(3,4)

      labels = [TPaveText( 0.6, 0.6, 0.9, 0.90, "NDC") for i in range(len(hNames)) ]

      for i, hName in enumerate(hNames):
        c.cd(i+1)
        hist = infile.Get(hName) 

        nSig, nSigError, fitStatus, frame = fitHist(hist)

        histDict[hName]["nSig"] = nSig
        histDict[hName]["nSigError"] = nSigError

        frame.Draw()

        labels[i].AddText("N Signal = " + str(nSig))
        labels[i].AddText("Err Sig  = " + str(nSigError))
        labels[i].AddText("Fit Status: " + str(fitStatus))
        labels[i].SetFillStyle(0)
        labels[i].Draw()

      c.SaveAs("fits/" + cName + ".png")
      c.Close()

for eVar in eVars:
  for eSel in eventSelections:
    h_tag = eSel + "invmass_" + eVar
    hNames = [h for h in histDict.keys() if h_tag in h]

    cName = dName + "_" + h_tag
    c = TCanvas(cName, cName, 900, 900)
      
    c.cd()
    c.Divide(3,4)

    labels = [TPaveText( 0.6, 0.6, 0.9, 0.9, 'NDC') for i in range(len(hNames))]

    for i, hName in enumerate(hNames):
      c.cd(i+1)
      hist = infile.Get(hName) 

      nSig, nSigError, fitStatus, frame = fitHist(hist)

      histDict[hName]["nSig"] = nSig
      histDict[hName]["nSigError"] = nSigError

      frame.Draw()

      labels[i].AddText("N Signal = " + str(nSig))
      labels[i].AddText("Err Sig  = " + str(nSigError))
      labels[i].AddText("Fit Status: " + str(fitStatus))
      labels[i].Draw()

    c.SaveAs("fits/" + cName + ".png")
    c.Close()


save_obj(histDict, dName+"_histDict")
