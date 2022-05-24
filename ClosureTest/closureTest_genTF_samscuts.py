from ROOT import TCanvas, TFile, TH1F, PyConfig
from ROOT import gROOT

#from utils import *
#from closureTestUtils import *
from collections import OrderedDict
import sys
import numpy as np


head = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/" 


dType = "WGJets"
era = "Summer16v3"
tag = "PtG"
HLT = 21

#dType = "TTJets"
#era = "Summer16v3"
#tag = "genMET-150"
#HLT = 21


versiontag = "samscuts"

outDir = "hists/"
###I/O
#ntuple_version =  "TreeMakerRandS_skimsv8"
ntuple_version = "TreeMaker"

ntuple_location = "root://cmsxrootd.fnal.gov//store/group/lpcsusyphotons/" + ntuple_version + "/"

#txtfiledir = "/scratch365/atownse2/SUSYDiPhoton/condor/" #I should make condor work here
txtfiledir = head + "Samples/"
txtfile = open(txtfiledir + ntuple_version + "_filelist.txt", "r")

filenames = [ ntuple_location + i[:-1] for i in txtfile if dType in i and tag in i and era in i]
###

### Some command line parsing to determine which files to run over
runbatch = len(sys.argv) > 1

if runbatch:
  batches = int(sys.argv[1])
  run     = int(sys.argv[2])

  file_batches = np.array_split(filenames, batches)

  filenames = file_batches[run]
  outDir = "/scratch365/atownse2/SUSYDiPhoton/ClosureTest/" + dType + "_" + era + "/"

if not any(filenames):
  print("No files in batch")
  quit()


print( "Running over " + str(len(filenames)) + " files")
###


###
outfilename = dType + "_" + era + "_" +  ntuple_version +  "_closureTest_genTF_ns" 

if versiontag:
  outfilename += "_" + versiontag
if runbatch:
  outfilename += "_" + str(run)

outfilename += ".root"

print( "outfilename : " +  outfilename)
print( "outfiledir  : " +  outDir)
###

genTransferFactor = {"barrel": (0.017383707625563526, 9.211294454590436e-05),
                     "endcap": (0.024732126508268247, 0.00017889773719544853)}

###


gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


hTypes = ["DiEMpt", "BDT", "met", "met_lowBDT", "met_medBDT", "met_highBDT", "leadpt", "trailpt", "leadeta", "traileta"]

selections = ["eg", "gg"]

fakerateRegions = ["barrel", "endcap"]

ptTags = ["tagHpt", "tagLpt"]

histDict = OrderedDict()

genTypes = ["genEle" , "all"]
recoTypes = ["recoPho", "recoEle"]

for gen in genTypes:
  for reco in recoTypes:
    hName = gen + "_" + reco
    histDict[hName] = TH1F(hName, hName, 1, 0, 1)



def getHistName(sel, hType, region = False, ptTag = False, bdt = False):
  hName = sel + "_" + hType + "_"
  if bdt:
    hName += bdt + "_"
  if region:
    hName += region + "_"
  if ptTag:
    hName += ptTag + "_"

  return hName

from closureTestUtils import varBins, deltaR

for hType in hTypes:
  for varName, b in varBins.items():
    if varName in hType:
      bins = b

  for sel in selections:
    if sel == "eg":
      for region in fakerateRegions:
        for tag in ptTags:
          hName = getHistName(sel, hType, region = region, ptTag = tag)
          histDict[hName] = TH1F(hName, hName, len(bins)-1, bins)
    else:
      hName = getHistName(sel, hType)
      histDict[hName] = TH1F(hName, hName, len(bins)-1, bins)



##Start Loop over files in batch
for filename in filenames:
 
  f = TFile.Open(filename)
  print("Opening File " + filename )
  t = f.Get('TreeMaker2/PreSelection')

  testCount = 0
  ##
  print("Starting Event Loop")
  for event in t:
    if testCount == 0:
      print("Made it into event loop")

    testCount += 1

    if t.TriggerPass[HLT]==1:

      #print("Passes HLT")

      #if t.IsRandS != 0:
      #  continue

      if len(t.Photons) < 2: 
        continue

      EMpass = [ {"pt"  : t.Photons[i].Pt(),
                  "eta" : t.Photons[i].Eta(),
                "4vec"  : t.Photons[i] , 
                "xseed" : bool(t.Photons_hasPixelSeed[i]),
                "barrel": t.Photons_isEB[i],
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 30 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

      if len(EMpass) < 2:
        continue

      em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first

      if not em[0]["barrel"]  and not em[1]["barrel"]:
        continue

      #jets = [j for j in t.hadronicJets if j.Pt() > 30]

      #nJets = len(jets)

      #vmult = t.NVtx

      MET = t.MET

      w = t.CrossSection*t.puWeight*1000*35.9

      #BDT_val = t.mva_BDT

      #BDT = getBDTbin(BDT_val)

      ### Just use two highest Pt objects
      lead = em[0]
      trail = em[1]      


      pcdiff = lambda a, b : 100*np.abs(a-b)/np.average([a,b])
      #Require that exactly one of thsese is matched to a gen Electron
      nmatch = 0
      for e in [lead, trail]:
  
        if e["xseed"]:
          histDict["all_recoEle"].Fill(0)
        else:
          histDict["all_recoPho"].Fill(0)

        for genEle in t.GenElectrons:
          if genEle.Pt() < 10:
            continue

          if deltaR(e["4vec"], genEle) < 0.2 and pcdiff(e["pt"], genEle.Pt()) < 20:
            nmatch += 1
            if e["xseed"]: #Is an electron
              histDict["genEle_recoEle"].Fill(0)
            else: #Is a photon
              histDict["genEle_recoPho"].Fill(0)

            break
      
      if nmatch != 1:
        continue

      if lead["xseed"] != trail["xseed"]: # e + gamma

        if lead["xseed"]: #Find the region/tag of the electron
          tag = "tagHpt"
          if lead["barrel"]:
            region = "barrel"
          else:
            region = "endcap"
        else:
          tag = "tagLpt"
          if trail["barrel"]:
            region = "barrel"
          else:
            region = "endcap"

        histDict[getHistName("eg", "met", region, tag)].Fill(MET, w)
        #histDict[getHistName("eg", "met", region, tag, bdt = BDT)].Fill(MET, w)

        histDict[getHistName("eg", "leadpt", region, tag)].Fill(lead["pt"], w) 
        histDict[getHistName("eg", "trailpt", region, tag)].Fill(trail["pt"], w)
        histDict[getHistName("eg", "leadeta", region, tag)].Fill(abs(lead["eta"]), w)
        histDict[getHistName("eg", "traileta", region, tag)].Fill(abs(trail["eta"]), w)

        DiEMpt = (lead["4vec"] + trail["4vec"]).M()

        histDict[getHistName("eg", "DiEMpt", region, tag)].Fill(DiEMpt, w)


      elif not (lead["xseed"] or trail["xseed"]): # gamma + gamma

        histDict[getHistName("gg", "met")].Fill(MET, w)
        #histDict[getHistName("gg", "met", bdt = BDT)].Fill(MET, w)

        histDict[getHistName("gg", "leadpt")].Fill(lead["pt"], w)
        histDict[getHistName("gg", "trailpt")].Fill(trail["pt"], w)
        histDict[getHistName("gg", "leadeta")].Fill(abs(lead["eta"]), w)
        histDict[getHistName("gg", "traileta")].Fill(abs(trail["eta"]), w)

        DiEMpt = (lead["4vec"] + trail["4vec"]).M()

        histDict[getHistName("gg", "DiEMpt")].Fill(DiEMpt, w)


      ###End of if TriggerPass
    ###End of event loop 

  f.Close()





print("Writing hists to: "+outDir+outfilename)
outFile = TFile.Open(outDir+outfilename, "RECREATE")
outFile.cd()



###Write Hists
for hName, hist in histDict.items():
  hist.Write()

outFile.Close()
