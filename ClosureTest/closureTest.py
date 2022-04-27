from ROOT import TCanvas, TFile, TH1F, PyConfig
from ROOT import gROOT

from utils import *
from closureTestUtils import *
from collections import OrderedDict

dType = "WGJets"
run = "Summer16v3"
ptTag = "PtG"
HLT = 21

#dType = "TTJets"
#run = "Summer16v3"
#ptTag = "genMET-150"
#HLT = 21


#txtfiledir = "/scratch365/atownse2/SUSYDiPhoton/condor/"
txtfiledir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/samples/"

#outDir = "/scratch365/atownse2/SUSYDiPhoton/closure/" + dType + "_" + run + "/"
outDir = "hists/"

outfilename = dType + "_" + run + "_closure"
txtfile = open(txtfiledir + "skimsv8_filelist.txt", "r")



print "outfilename : " +  outfilename
print "outfiledir  : " +  outDir

filenames = ["root://cmsxrootd.fnal.gov//store/group/lpcsusyphotons/TreeMakerRandS_skimsv8/"+ i[:-1] for i in txtfile if dType in i and run in i and ptTag in i]
#filenames = ["root://cmsxrootd.fnal.gov//store/group/lpcsusyphotons/TreeMaker/"+ i[:-1] for i in txtfile if dType in i and "Summer16v3" in i]

print(filenames[0])

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


hTypes = ["DiEMpt", "BDT", "met", "met_lowBDT", "met_medBDT", "met_highBDT", "leadpt", "trailpt", "leadeta", "traileta"]

# DiEMpt  = magnitude of vector sum of pt

selections = ["eg", "gg"]

fakerateRegions = ["barrel", "endcap"]

ptTags = ["tagHpt", "tagLpt"]

histDict = OrderedDict()

ptRanges = ptRangeDict[dType]

for ptRange in ptRanges:

  hName = "hEvents_" + ptRange
  histDict[hName] = TH1F(hName, hName, 1, 0, 1 )

  for hType in hTypes:
    for varName, b in varBins.items():
      if varName in hType:
        bins = b

    for sel in selections:
      if sel == "eg":
        for region in fakerateRegions:
          for tag in ptTags:
            hName = getHistName(sel, hType, ptRange, region = region, ptTag = tag)
            histDict[hName] = TH1F(hName, hName, len(bins)-1, bins)
      else:
        hName = getHistName(sel, hType, ptRange)
        histDict[hName] = TH1F(hName, hName, len(bins)-1, bins)


##Start Loop over files in batch
for filename in filenames[nStart:nEnd]:
 
  f = TFile.Open(filename)
  print("Opening File " + filename )
  t = f.Get('TreeMaker2/PreSelection')

  for ptRange in ptRanges:
    if ptRange in filename:
      ptTag = ptRange
      break

  if not ptTag:
    print("Tag not found:" + filename)
    quit()

  tCounter = f.Get("tcounter")

  for i in range(tCounter.GetEntries()):
    histDict["hEvents_"+ptTag].Fill(1)

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
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 80 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

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

      BDT_val = t.mva_BDT

      BDT = getBDTbin(BDT_val)

      ### Just use two highest Pt objects
      lead = em[0]
      trail = em[1]      

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

        histDict[getHistName("eg", "met", ptTag, region, tag)].Fill(MET, w)
        histDict[getHistName("eg", "met", ptTag, region, tag, bdt = BDT)].Fill(MET, w)

        histDict[getHistName("eg", "leadpt", ptTag, region, tag)].Fill(lead["pt"], w) 
        histDict[getHistName("eg", "trailpt", ptTag, region, tag)].Fill(trail["pt"], w)
        histDict[getHistName("eg", "leadeta", ptTag, region, tag)].Fill(abs(lead["eta"]), w)
        histDict[getHistName("eg", "traileta", ptTag, region, tag)].Fill(abs(trail["eta"]), w)

        DiEMpt = (lead["4vec"] + trail["4vec"]).M()

        histDict[getHistName("eg", "DiEMpt", ptTag, region, tag)].Fill(DiEMpt, w)


      elif not (lead["xseed"] or trail["xseed"]): # gamma + gamma

        histDict[getHistName("gg", "met", ptTag)].Fill(MET, w)
        histDict[getHistName("gg", "met", ptTag, bdt = BDT)].Fill(MET, w)

        histDict[getHistName("gg", "leadpt", ptTag)].Fill(lead["pt"], w)
        histDict[getHistName("gg", "trailpt", ptTag)].Fill(trail["pt"], w)
        histDict[getHistName("gg", "leadeta", ptTag)].Fill(abs(lead["eta"]), w)
        histDict[getHistName("gg", "traileta", ptTag)].Fill(abs(trail["eta"]), w)

        DiEMpt = (lead["4vec"] + trail["4vec"]).M()

        histDict[getHistName("gg", "DiEMpt", ptTag)].Fill(DiEMpt, w)


      ###End of if TriggerPass
    ###End of event loop 

  f.Close()


foutName = outDir + outfilename + "_" + str(nBatchRun) + ".root"

print(foutName)

outFile = TFile.Open(foutName, "RECREATE")
outFile.cd()



###Write Hists
for hName, hist in histDict.items():
  hist.Write()

outFile.Close()
