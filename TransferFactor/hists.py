from collections import OrderedDict
from ROOT import TFile, TH1F
import numpy as np
from utils import l2b

def get_hist_name(selection, region, var = False, var_bin = False, tag = False,  bdt = False, ptRange = False):
  hName = selection + "invmass_" + region

  if var and var_bin:
    hName += "_" + var + str(var_bin[0]) + "to" + str(var_bin[1]) 
  if tag:
    hName += "_" + tag
  if var:
    hName += "_" + var
  if bdt:
    hName += "_" + bdt
  if ptRange:
    hName += "_" + ptRange

  return hName

#Define and initialize hists
ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets_DiLept_TuneCUETP8M1": ["_"]}

varBins = { "pt"    : [80,100,150,200,300],
            "eta"   : np.linspace(0,2.4, num = 8),
            "jnum"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250],
             "all"   : ["all"] }

tags = ["tagHpt", "tagLpt", "all"]

BDTbins = ["lowBDT", "medBDT", "highBDT"]

kVars = ["pt", "eta"]

eVars = ["jnum", "vmult", "met"]

regions = ["BB" , "BE", "EB", "EE"]

eventSelections = ["ee","eg"]

histDict = OrderedDict()

def init_hists(dType, useBDTbins = False, init_scaled_eg = False, isMC = False):

  for region in regions:
    for sel in eventSelections:
      hName = get_hist_name( sel, region)
      h = TH1F(hName, hName, 128, 0, 150)
      h.SetDirectory(0)
      histDict[hName] = h

  """
  for var in kVars+eVars:
    for var_bin in l2b(varBins[var]):
      for region in regions:
        for sel in eventSelections:
          hName = get_hist_name(sel, region, var, var_bin)
          h = TH1F(hName, hName, 128, 0, 150)
          h.SetDirectory(0)
          histDict[hName] = h
  """
  ### Add new hists here
  

  ###

  if isMC:
    ptRanges = ptRangeDict[dType]

    for hName, hist in histDict.items():
      for ptRange in ptRanges:
        hName_ptRange = hName + "_" + ptRange
        h = TH1F(hName_ptRange, hName_ptRange, 128, 0, 150)
        h.SetDirectory(0)
        histDict[hName_ptRange] = h

      del histDict[hName]

    for ptRange in ptRanges:
      hName = "hEvents_" + ptRange
      h = TH1F(hName, hName, 1, 0, 1)
      histDict[hName] = h

  return


def fill_hists(sel, region, invmass, ptRange=False, **varDict):

  histDict[get_hist_name(sel, region, ptRange=ptRange)].Fill(invmass)
  
  """
  for var, value in varDict.items(): ##For this to work I need to seperate the objects by Pt which involves many changes
    #Determine which bin
    bins = l2b(varBins[var])
    foundBin = False #TODO include overflows?
    for var_bin in bins:
      if value > var_bin[0] and value <= var_bin[1]:
        foundBin = True
        break

    if foundBin:
      histDict[get_hist_name(sel, region, var, var_bin)].Fill(invmass)
  """

  return

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
    h = TH1F(hName, hName, 128, 0, 150)
    h.SetDirectory(0)

    for ptRange in ptRanges:
      h_events = histDict["hEvents_" + ptRange]
      h_ptRange = histDict[hName+"_"+ptRange]

      h_ptRange.Scale(1000000/h_events.GetEntries())
      h.Add(h_ptRange)

    scaled_hists[hName] = h

  return

def write_hists(fout_name):
  outFile = TFile.Open(fout_name, "RECREATE")
  outFile.cd()

  ###Write Hists
  for hName, hist in histDict.items():
    hist.Write()

  outFile.Close()
  return
