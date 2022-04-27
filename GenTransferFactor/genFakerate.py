from ROOT import TCanvas, TFile, TH1F, PyConfig
from ROOT import gROOT
import sys
import numpy as np
from utils import dR
from array import array
from collections import OrderedDict

#dType = "DYJetsToLL_M-50"
#run = "Summer16v3"
#rtag = "HT"
#HLT = 21

dType = "WGJets"
run = "Summer16v3"
rtag = "PtG"
HLT = 21

versiontag = ""

#txtfiledir = "/scratch365/atownse2/SUSYDiPhoton/condor/"
txtfiledir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/samples/"

#outDir = "/scratch365/atownse2/SUSYDiPhoton/fakerate/" + dType + "_" + run + "/"
outDir = "hists/"

outfilename = dType + "_" + run + "_genfakerate_ns" 

if versiontag:
  outfilename += "_" + versiontag

txtfile = open(txtfiledir + "skimsforAustin_filelist.txt", "r")

print "outfilename : " +  outfilename
print "outfiledir  : " +  outDir


filenames = ["root://cmsxrootd.fnal.gov//store/group/lpcsusyphotons/skimsforAustin/"+ i[:-1] for i in txtfile if dType in i and rtag in i and run in i]


print("Number of files to run over: " + str(len(filenames)))

### Some command line parsing to determine which files to run over
nBatch = int(sys.argv[1])
nBatchRun = int(sys.argv[2])

nSize = len(filenames)//nBatch + 1
nStart = nSize*nBatchRun
nEnd = nSize*(nBatchRun + 1)

if nEnd > len(filenames):
  nEnd = len(filenames)

if nStart > len(filenames):
  print("No more batches")
  quit()

print("nStart = " + str(nStart))
print("nEnd = " + str(nEnd))

###


gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


### Define hist dicts and bins
ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets_DiLept_TuneCUETP8M1": ["_"]}

varBins = { "pt"    : [80,100,120,150,200,300],
            "eta"   : np.linspace(-2.4, 2.4, num=24),
            "njets"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250]}

genTypes = ["genEle"]

recoTypes = ["recoPho", "recoEle"]

regions = ["barrel" , "endcap"]

eVars = ["pt", "eta" , "njets" , "met"]

ptRanges = ptRangeDict[dType]

def getHistName( ptRange, genType, recoType, region, var):
  return genType + "_" + recoType + "_" + region + "_" + var + "_" + ptRange

hists = OrderedDict()

for ptRange in ptRanges:
  hists["hEvents_"+ptRange] = TH1F("hEvents_"+ptRange, "hEvents_"+ptRange, 1, 0, 1 ) 

  for genType in genTypes:
    for recoType in recoTypes:
      for region in regions:
        for eVar in eVars:
          h_name = getHistName(ptRange, genType, recoType, region, eVar)
          h_bins = array("d", varBins[eVar])
          hists[h_name] = TH1F( h_name, h_name, len(h_bins)-1, h_bins)

###

### Event counters (currently for testing only)
eventCount = 0
triggerPassCount = 0
triggerPass2EMCount = 0


##Start Loop over files in batch
for filename in filenames[nStart:nEnd]:
  print(filename)
 
  f = TFile.Open(filename)
  print("Opening File " + filename )
  t = f.Get('TreeMaker2/PreSelection')

  #print(t.ls())

  for ptRange in ptRanges:
    if ptRange in filename:
      ptTag = ptRange
      break

  tCounter = f.Get("tcounter")

  for i in range(tCounter.GetEntries()):
    hists["hEvents_"+ptTag].Fill(1)

  ##
  print("Starting Event Loop")
  for event in t:
    eventCount += 1

    if t.TriggerPass[HLT]==1:

      #print("Passes HLT")

      triggerPassCount += 1

      if len(t.Photons) < 2: 
        continue

      EMpass = [ {"pt"  : t.Photons[i].Pt(),
                "4vec"  : t.Photons[i] , 
                "xseed" : bool(t.Photons_hasPixelSeed[i]),
                "barrel": t.Photons_isEB[i],
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 80 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

      if len(EMpass) < 2:
        continue

      em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first

      if not t.Photons_isEB[em[0]["index"]]  and not t.Photons_isEB[em[1]["index"]]:
        continue

      jets = [j for j in t.hadronicJets if j.Pt() > 30]

      """
      for jet in [j for j in t.Jets if j.Pt() > 30]:
        jetPass = True
        for emObj in em:
          if dR(jet,emObj[1]) < 0.4:
            jetPass = False
        if jetPass:
          jets.append(jet)
      """

      njets = len(jets)

      vmult = t.NVtx

      met = t.MET

      w = t.CrossSection*t.puWeight*1000*35.9

      for reco in em:
        if reco["barrel"]:
          region = "barrel"
        else:
          region = "endcap"

        for genEle in t.GenElectrons:
          if genEle.Pt() < 10:
            continue

          pcdiff = lambda a, b : 100*np.abs(a-b)/np.average([a,b])

          if dR(reco["4vec"], genEle) < 0.03 and pcdiff(reco["pt"], genEle.Pt()) < 20:
            
            if reco["xseed"]: #Electron
              #Fill Hists
              hists[getHistName(ptRange, "genEle", "recoEle", region, "pt" )].Fill(reco["4vec"].Pt() , w)
              hists[getHistName(ptRange, "genEle", "recoEle", region, "eta")].Fill(reco["4vec"].Eta(), w)
              hists[getHistName(ptRange, "genEle", "recoEle", region, "njets")].Fill(njets, w)
              hists[getHistName(ptRange, "genEle", "recoEle", region, "met")].Fill(met, w)

            else: #Electron faking photon
              #Fill Hists
              hists[getHistName(ptRange, "genEle", "recoPho", region, "pt" )].Fill(reco["4vec"].Pt() , w)
              hists[getHistName(ptRange, "genEle", "recoPho", region, "eta")].Fill(reco["4vec"].Eta(), w)
              hists[getHistName(ptRange, "genEle", "recoPho", region, "njets")].Fill(njets, w)
              hists[getHistName(ptRange, "genEle", "recoPho", region, "met")].Fill(met, w)

            break
      
      ###End of if TriggerPass
    ###End of event loop 
  f.Close()

print("Event count = " + str(eventCount))
print("Trigger Pass count = " + str(triggerPassCount))
print("Trigger Pass count w/ 2 em objects = " + str(triggerPass2EMCount))



foutName = outDir + outfilename + "_" + str(nBatchRun) + ".root"

print(foutName)

outFile = TFile.Open(foutName, "RECREATE")
outFile.cd()

###Write Hists
for h_name, hist in hists.items():
  hist.Write()


