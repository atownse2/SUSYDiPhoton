from collections import OrderedDict
from ROOT import TFile, TH1F
from array import array

def l2b(l):
  b = []
  for i in range(len(l)-1):
    b.append((l[i],l[i+1]))
  return b

def getHistName(eventSelection, var_name, ptRange, tag = None, BDTbin = None, scaled_eg = False):

  if scaled_eg:
    hName = "scaled_" + eventSelection + "_"
  else:
    hName = eventSelection + "_" 
  if tag:
    hName += tag + "_"

  if BDTbin:
    hName += BDTbin + "_"

  return hName + var_name + "_" + ptRange


ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets_DiLept_TuneCUETP8M1": ["_"]}

varBins = { "pt"    : [80,100,150,200,300],
            "eta"   : ["barrel", "endcap"],
            "jnum"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,10,20,35,50,70,90,120,200],
             "all"   : ["all"] }

eventSelections = ["eg", "gg"]

BDTbins = ["lowBDT", "medBDT", "highBDT"]

eVars = [ "met"]

histDict = OrderedDict()

def init_hists(dType):

  ptRanges = ptRangeDict[dType]

  for ptRange in ptRanges:
    hName = "hEvents_" + ptRange
    histDict[hName] = {}
    histDict[hName]["hist"] = TH1F(hName, hName, 1, 0, 1 )


  for eSel in eventSelections:
    for eVar in eVars:
      bins = array("d", varBins[eVar])
      for ptRange in ptRanges:
        hName = getHistName(eSel, eVar, ptRange)
        histDict[hName] = {}
        histDict[hName]["hist"] = TH1F(hName, hName, len(bins)-1, bins)
        for BDTbin in BDTbins:
          hName = getHistName(eSel, eVar, ptRange, BDTbin = BDTbin)
          histDict[hName] = {}
          histDict[hName]["hist"] = TH1F(hName, hName, len(bins)-1, bins)

  for hName in [h for h in histDict.keys() if "eg_" in h]:
    scaledName = hName.replace("eg_", "scaled_eg_")
    histDict[scaledName] = {}
    histDict[scaledName]["hist"] = histDict[hName]["hist"].Clone(scaledName)

  return


def fill_hist(ptRange, eventSelection, var_name, var, weight, tag = None, BDTbin = None, scaled_eg = False):
  hName = getHistName(eventSelection, var_name, ptRange, tag = tag, BDTbin = BDTbin, scaled_eg = scaled_eg)
  histDict[hName]["hist"].Fill(var, weight)
  return
