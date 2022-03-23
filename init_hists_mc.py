from collections import OrderedDict
from ROOT import TFile, TH1F

def l2b(l):
  b = []
  for i in range(len(l)-1):
    b.append((l[i],l[i+1]))
  return b

def getHistName(eventSelection, var_name, var_bin, ptRange, tag = None, BDTbin = None, scaled_eg = False):

  if scaled_eg:
    hName = "scaled_" + eventSelection + "invmass_" + ptRange + "_"
  else:
    hName = eventSelection + "invmass_" + ptRange + "_"

  if tag:
    hName += tag + "_"

  if BDTbin:
    hName += BDTbin + "_"

  if "eta" in var_name:
    return hName + var_name + "_" + var_bin
  else:
    return hName + var_name + "_" + str(var_bin[0]) + "to" + str(var_bin[1]) 


useBDTbins = False
use_scaled_eg = False

ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets_DiLept_TuneCUETP8M1": ["_"]}

varBins = { "pt"    : [80,100,150,200,300],
            "eta"   : ["barrel", "endcap"],
            "jnum"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250],
             "all"   : ["all"] }

eventSelections = ["ee", "eg", "gg"]

tags = ["tagHpt", "tagLpt", "all"]

BDTbins = ["lowBDT", "medBDT", "highBDT"]

kVars = ["pt", "eta"]

eVars = ["jnum", "vmult", "met"]

t_hists = ["eeinvmass","eginvmass", "gginvmass"]


histDict = OrderedDict()

def init_hists(dType, useBDTbins = False, init_scaled_eg = False):

  ptRanges = ptRangeDict[dType]

  for ptRange in ptRanges:
    hName = "hEvents_" + ptRange
    histDict[hName] = {}

  for t in t_hists:
    for ptRange in ptRanges:
      hName = t + "_" + ptRange   
      histDict[hName] = {}

  for eSel in eventSelections:
    for kVar in kVars:
      for ptRange in ptRanges:

        if "eta" in kVar:
          bins = varBins["eta"]
        else:
          bins = l2b(varBins[kVar])
        for b in bins:
          for tag in tags:
            hName = getHistName(eSel, kVar, b, ptRange, tag = tag)
            histDict[hName] = {}


  for eSel in eventSelections:
    for eVar in eVars:
      for ptRange in ptRanges:
        bins = l2b(varBins[eVar])
        for b in bins:
          if eVar == "met" and useBDTbins:
            hName = getHistName(eSel, eVar, b, ptRange)
            histDict[hName] = {}
            for BDTbin in BDTbins:
              hName = getHistName(eSel, eVar, b, ptRange, BDTbin = BDTbin)
              histDict[hName] = {}
          else:
            hName = getHistName(eSel, eVar, b, ptRange)
            histDict[hName] = {}

  if init_scaled_eg:
    use_scaled_eg = True
    for hName in [h for h in histDict.keys() if "eginvmass" in h]:
      scaledName = hName.replace("eginvmass", "scaled_eginvmass")
      histDict[scaledName] = {}


  for hName, hDict in histDict.items():

    if "hEvents" in hName:
      hDict["hist"] = TH1F(hName, hName, 1, 0, 1 )
    else:
      hDict["hist"] = TH1F(hName, hName, 128, 0, 150)


def fill_hist(ptRange, eventSelection, var_name, var, invMass, weight, tag = None, BDTbin = None, scaled_eg = False):
  var_bin = False

  if "eta" in var_name:
    if var:
      var_bin = "barrel"
    else:
      var_bin = "endcap"
  else:
    bins = l2b(varBins[var_name])
    for b in bins:
      if var >= b[0] and var < b[1]:
        var_bin = b
        break
    if not var_bin and var >= bins[-1][1]: ##Includes overflow in last bin
      var_bin = bins[-1]
    if not var_bin and var < bins[0][0]:
      var_bin = bins[0]
  
  if var_bin:
    hName = getHistName(eventSelection, var_name, var_bin, ptRange, tag = tag, BDTbin = BDTbin)
    histDict[hName]["hist"].Fill(invMass, weight)

  else:
    print("Not correctly setting var bin:" + var_name + " " + str(var))
    quit()
