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


hTypes = ["met", "met_lowBDT", "met_medBDT", "met_highBDT", "leadpt", "trailpt", "leadeta", "traileta"]

selections = ["scaled_eg", "gg"]

fakerateRegions = ["barrel", "endcap"]

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
      if sel == "scaled_eg":
        for region in fakerateRegions:
          hName = getHistName(sel, hType, ptRange, region)
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
                "xseed" : t.Photons_hasPixelSeed[i],
                "barrel": t.Photons_isEB[i],
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 80 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

      if len(EMpass) != 2: #TODO change here
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

      phoPass = [i for i in em if not i["xseed"]]

      elePass = [i for i in em if i["xseed"]]

      ### Just compare selections with exactly 2 EM objects

      if len(phoPass) >= 1 and len(elePass) >= 1:

        for ele in elePass:
          for pho in phoPass:
            if dR(ele["4vec"],pho["4vec"]) < 0.4:
              continue

            if ele["barrel"]:
              region = "barrel"
            else:
              region = "endcap"
 
            histDict[getHistName("scaled_eg", "met", ptTag, region)].Fill(MET, w)
            histDict[getHistName("scaled_eg", "met", ptTag, region, bdt = BDT)].Fill(MET, w)

            if ele["pt"] > pho["pt"]:
              lead = ele
              trail = pho
            else:
              lead = pho
              trail = ele

            histDict[getHistName("scaled_eg", "leadpt", ptTag, region)].Fill(lead["pt"], w) #TODO check combinatorics
            histDict[getHistName("scaled_eg", "trailpt", ptTag, region)].Fill(trail["pt"], w)
            histDict[getHistName("scaled_eg", "leadeta", ptTag, region)].Fill(lead["eta"], w)
            histDict[getHistName("scaled_eg", "traileta", ptTag, region)].Fill(trail["eta"], w)


      if len(phoPass) >= 2: 

        lead = phoPass[0]
        trail = phoPass[1]

        histDict[getHistName("gg", "met", ptTag)].Fill(MET, w)
        histDict[getHistName("gg", "met", ptTag, bdt = BDT)].Fill(MET, w)


        histDict[getHistName("gg", "leadpt", ptTag)].Fill(lead["pt"], w)
        histDict[getHistName("gg", "trailpt", ptTag)].Fill(trail["pt"], w)
        histDict[getHistName("gg", "leadeta", ptTag)].Fill(lead["eta"], w)
        histDict[getHistName("gg", "traileta", ptTag)].Fill(trail["eta"], w)


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
