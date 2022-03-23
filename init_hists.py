from collections import OrderedDict
from ROOT import TFile, TH1F

def l2b(l):
  b = []
  for i in range(len(l)-1):
    b.append((l[i],l[i+1]))
  return b

def getHistName(eventSelection, var_name, var_bin, tag):
  hName = eventSelection + "invmass_"

  if tag:
    hName += tag + "_"

  if "eta" in var_name:
    return hName + var_name + "_" + var_bin
  else:
    return hName + var_name + "_" + str(var_bin[0]) + "to" + str(var_bin[1]) 


varBins = { "pt"    : [80,100,150,200,300],
            "eta"   : ["barrel", "endcap"],
            "jnum"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250],
             "all"   : ["all"] }


eventSelections = ["ee", "eg"]

tags = ["tagHpt", "tagLpt", "all"]

kVars = ["pt", "eta"]

eVars = ["jnum", "vmult", "met"]



histDict = OrderedDict({"eeinvmass": {},
         "eginvmass": {},
         "gginvmass": {}})

for eSel in eventSelections:
  for kVar in kVars:
    if "eta" in kVar:
      bins = varBins["eta"]
    else:
      bins = l2b(varBins[kVar])
    for b in bins:
      for tag in tags:
        hName = getHistName(eSel, kVar, b, tag)
        histDict[hName] = {}

for eSel in eventSelections:
  for eVar in eVars:
    bins = l2b(varBins[eVar])
    for b in bins:
      hName = getHistName(eSel, eVar, b, None)
      histDict[hName] = {}



def init_hists():
  print("Initializing histograms")
  for hName, hDict in histDict.items():
    hDict["hist"] = TH1F(hName, hName, 128, 0, 150)


def fill_hist(eventSelection, tag, var_name, var, invMass):
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
  
  if var_bin:
    hName = getHistName(eventSelection, var_name, var_bin, tag)
    histDict[hName]["hist"].Fill(invMass)
  else:
    print("Not correctly setting var bin:" + var_name + " " + str(var))
    quit()



