from array import array
from ROOT import TFile, TH1F
from collections import OrderedDict
import numpy as np

def getHistName(sel, hType, region = False, ptTag = False, bdt = False, ptRange = False):
  hName = sel + "_" + hType + "_"

  if bdt:
    hName += bdt + "_"

  if region:
    hName += region + "_"

  if ptTag:
    hName += ptTag + "_"

  if ptRange:
    hName += ptRange

  return hName


def cloneDimsTH1F(hName, hist):
  tmp = TH1F(hName, hName, hist.GetNbinsX(), hist.GetXaxis().GetXbins().GetArray())
  return tmp



### Define hist dicts and bins
ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets": ["SingleLeptFromTbar", "SingleLeptFromT", "DiLept"]}

varBins = { "pt"    : array("d", np.round(np.linspace(80,300, num = 32), decimals = 3)),
            "eta"   : array("d", np.round(np.linspace(0,2.4, num = 32), decimals = 2)),
            "met"   : array("d", [0,10,20,35,50,70,90,110,130,150,150,185,185,250]),
            "DiEMpt": array("d", np.round(np.linspace(0, 500, num = 64), decimals = 1)),
            "BDT"   : array("d", np.round(np.linspace(0, 1  , num = 64), decimals = 1))}

def get_bins(hType):
  return varBins[hType.split("_")[0]]

### Initialize and load hists

hTypes = ["DiEMpt", "BDT", "met", "met_lowBDT", "met_medBDT", "met_highBDT", "pt_lead", "pt_trail", "eta_lead", "eta_trail"]

# DiEMpt  = magnitude of vector sum of pt

selections = ["eg", "gg"]

regions = ["barrel", "endcap"]

ptTags = ["tagHpt", "tagLpt"]

histDict = OrderedDict()


### Initialize histograms
def init_hists(dType):
  ptRanges = ptRangeDict[dType]
  for ptRange in ptRanges:
    
    hName = "hEvents_" + ptRange
    histDict[hName] = TH1F(hName, hName, 1, 0, 1 )

    for hType in hTypes:
      bins = get_bins(hType)

      for sel in selections:
        if sel == "eg":
          for region in regions:
            for tag in ptTags:
              hName = getHistName(sel, hType, region = region, ptTag = tag, ptRange = ptRange)
              histDict[hName] = TH1F(hName, hName, len(bins)-1, bins)
        else:
          hName = getHistName(sel, hType, ptRange = ptRange)
          histDict[hName] = TH1F(hName, hName, len(bins)-1, bins)

  for hName, hist in histDict.items():
    hist.SetDirectory(0) #Write to RAM

  return histDict


def load_hists(dType, filename):
  f = TFile.Open(filename)

  if not f:
    print("File unable to be opened: " + filename)
    quit()

  for key in f.GetListOfKeys():
    h = key.ReadObj()
    h.SetDirectory(0) #Write to RAM
    histDict[h.GetName()] = h

  f.Close()
  return

scaled_hists = OrderedDict()

def scale_hists(dType):
  ##Find list of unique hist names:
  ptRanges = ptRangeDict[dType]

  hNames = []

  for hName in histDict.keys():
    if "hEvents" in hName:
      continue
    if ptRanges[0] in hName:
      hNames.append(hName.replace("_"+ptRanges[0],""))

  for hName in hNames:
    h = cloneDimsTH1F(hName, histDict[hName+"_"+ptRanges[0]])
    h.SetDirectory(0)

    for ptRange in ptRanges:
      h_events = histDict["hEvents_" + ptRange]
      h_ptRange = histDict[hName+"_"+ptRange]

      h_ptRange.Scale(1000000/h_events.GetEntries())
      h.Add(h_ptRange)

    scaled_hists[hName] = h

  return

def write_hists(fout_name):
  print("Writing hists to: " + fout_name)

  outFile = TFile.Open(fout_name, "RECREATE")
  outFile.cd()

  ###Write Hists
  for hName, hist in histDict.items():
    hist.Write()

  outFile.Close()
  return
