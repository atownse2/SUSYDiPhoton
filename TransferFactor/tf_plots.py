from ROOT import *
from utils import *
import os, sys, math
import numpy as np

from hists import load_hists, scale_hists, eventSelections, regions, get_hist_name, scaled_hists 

gStyle.SetOptStat(0)
gROOT.SetBatch(1)

dType = "DYJetsToLL_M-50"
run = "Summer16v3"

#dType = "DoubleEG"
#run = "2016"

#dType = "DoubleEG"
#run  = "2017"

#dType = "EGamma"
#run = "2018"

###I/O
script_name = "tf_MC"

version = "071822v1"

#ntuple_version =  "TreeMakerRandS_skimsv8"
ntuple_version = "skimsforAustin"

fDir = "hists/"

fName = get_file_name(dType, run, ntuple_version, script_name, version)

f = TFile(fDir + fName)

##Initialize and load hists
load_hists(dType, fDir+fName)
scale_hists(dType)


fitDict = {}

for hName, hist in scaled_hists.items():

  print(hName + " nEntries : " + str(hist.GetEntries()))

  fitDict[hName] = {}

  cName = dType + "_Run" + run + "_" + hName
  c = TCanvas("c", cName, 600, 600) 
  c.cd()

  nSig, nSigError, fitStatus, frame = fitHist(hist)

  fitDict[hName]["nSig"] = nSig
  fitDict[hName]["nSigError"] = nSigError

  #Draw and save
  frame.Draw()

  label = TPaveText( 0.6, 0.3, 0.9, 0.80, 'NDC') 
  label.AddText("N Signal = " + str(nSig))
  label.AddText("Err Sig  = " + str(nSigError))
  label.AddText("Fit Status: " + str(fitStatus))
  label.Draw()

  c.SaveAs("fits/" + cName + ".png")
  c.Close()

  print(" ")

"""
##Calculate the transfer factor
for region in regions:

  eeDict = fitDict[get_hist_name("ee", region)]
  egDict = fitDict[get_hist_name("eg", region)]

  eeSig = eeDict["nSig"]
  eeSigError = eeDict["nSigError"]

  egSig = egDict["nSig"]
  egSigError = egDict["nSigError"]

  tf, e = transferFactor(eeSig, eeSigError, egSig, egSigError)

  print(region + " Transfer Factor: " + str(np.round(tf, decimals=4)) + " error: " + str(np.round(e, decimals=4)))
"""

regions1 = ["BB", "EE", "EB", "BE"]
##Calculate the transfer factor
tfDict = {run : { "barrel" : {"tagHpt" : tuple() ,
                              "tagLpt" : tuple()},
                  "endcap" : {"tagHpt" : tuple(),
                              "tagLpt" : tuple()} }}

for region in regions1:
  print("Finding transfer factor in region: " + region)

  eeDict = fitDict[get_hist_name("ee", region)]
  egDict = fitDict[get_hist_name("eg", region)]

  eeSig = eeDict["nSig"]
  eeSigError = eeDict["nSigError"]

  egSig = egDict["nSig"]
  egSigError = egDict["nSigError"]

  tf, e = transferFactor(eeSig, eeSigError, egSig, egSigError)

  print( "region: " + region + " tf : " +str(tf))

  tf_region = "barrel" if region[1] == "B" else "endcap"

  fac = 1


  tfDict[run][tf_region]["tagHpt"] = (tf/fac,e/fac)
  tfDict[run][tf_region]["tagLpt"] = (tf/fac,e/fac)


print(tfDict)
