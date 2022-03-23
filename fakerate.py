from ROOT import TCanvas, TFile, TH1F, PyConfig
from ROOT import gROOT
import sys
import numpy as np
import utils

dType = "DoubleEG"
run = "2016"
HLT = 21

#dType = "DoubleEG"
#run  = "2017"
#HLT = 23

#dType = "EGamma"
#run = "2018"
#HLT = 22

#txtfiledir = "/scratch365/atownse2/SUSYDiPhoton/condor/"
txtfiledir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/samples/"

#outDir = "/scratch365/atownse2/SUSYDiPhoton/fakerate/" + dType + "_" + run + "/"
outDir = "hists/"

outfilename = dType + "_" + run + "_fakerate"
txtfile = open(txtfiledir + "skimsforAustin_filelist.txt", "r")

print "outfilename : " +  outfilename
print "outfiledir  : " +  outDir

filenames = ["root://cmsxrootd.fnal.gov//store/group/lpcsusyphotons/skimsforAustin/"+ i[:-1] for i in txtfile if (dType in i and run in i)]

#print(filenames)

print("Number of files to run over: " + str(len(filenames)))

###Parse Command Line arguments to determine batch size (nBatches, batchNum)
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



###

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


#Create hist dict structure and initialize hists
from init_hists import *
init_hists()


### Event counters (currently for testing only)
eventCount = 0
triggerPassCount = 0
triggerPass2EMCount = 0

fCounter = 0

##Start Loop over files in batch
for filename in filenames[nStart:nEnd]:
 
  fCounter += 1

  f = TFile.Open(filename)
  print("Opening File " + filename )

  #print("Opening file " + str(fCounter) + " out of " + str(fNum))

  t = f.Get('TreeMaker2/PreSelection')

  #t.Print()

  #print("Starting Event Loop")
  for event in t:

    eventCount += 1

    if t.TriggerPass[HLT]==1:

      #print("Passes HLT")

      triggerPassCount += 1

      if len(t.Photons) < 2: 
        continue

      EMpass = [ {"pt"  : t.Photons[i].Pt(),
                "4vec"  : t.Photons[i] , 
                "xseed" : t.Photons_hasPixelSeed[i],
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

      nJets = len(jets)

      vmult = t.NVtx

      MET = t.MET

      phoPass = [i for i in em if not i["xseed"]]

      elePass = [i for i in em if i["xseed"]]

      ###

      if len(elePass) >= 2:
        if dR(elePass[0]["4vec"], elePass[1]["4vec"]) < 0.4:
          continue

        invMass = (elePass[0]["4vec"] + elePass[1]["4vec"]).M()
        histDict["eeinvmass"]["hist"].Fill(invMass)

        fill_hist("ee", "tagHpt", "pt" , elePass[1]["pt"]    , invMass)
        fill_hist("ee", "tagHpt", "eta", elePass[1]["barrel"], invMass)

        fill_hist("ee", "tagLpt", "pt" , elePass[0]["pt"]    , invMass)       
        fill_hist("ee", "tagLpt", "eta", elePass[0]["barrel"], invMass)

        fill_hist("ee", "all", "pt" , elePass[0]["pt"]    , invMass)       
        fill_hist("ee", "all", "eta", elePass[0]["barrel"], invMass)

        fill_hist("ee", "all", "pt" , elePass[1]["pt"]    , invMass)       
        fill_hist("ee", "all", "eta", elePass[1]["barrel"], invMass)

        fill_hist("ee" , None, "jnum" , nJets, invMass )
        fill_hist("ee" , None, "vmult", vmult, invMass )
        fill_hist("ee" , None, "met"  , MET  , invMass )

      
      if len(phoPass) >= 1 and len(elePass) >= 1:
        for ele in elePass:
          for pho in phoPass:
            if dR(ele["4vec"],pho["4vec"]) < 0.4:
              continue
            invMass = (ele["4vec"] + pho["4vec"]).M()
            histDict["eginvmass"]["hist"].Fill(invMass)

            if ele["pt"] > pho["pt"]:
              fill_hist("eg", "tagHpt", "pt" , pho["pt"]    , invMass)
              fill_hist("eg", "tagHpt", "eta", pho["barrel"], invMass)
            else:
              fill_hist("eg", "tagLpt", "pt" , pho["pt"]    , invMass)   
              fill_hist("eg", "tagLpt", "eta", pho["barrel"], invMass)

            fill_hist("eg", "all", "pt" , pho["pt"]    , invMass)   
            fill_hist("eg", "all", "eta", pho["barrel"], invMass)

            fill_hist("eg", None, "jnum" , nJets, invMass)
            fill_hist("eg", None, "vmult", vmult, invMass)
            fill_hist("eg", None, "met"  , MET  , invMass)

      if len(phoPass) >= 2:   
        invMass = (phoPass[0]["4vec"] + phoPass[1]["4vec"]).M()
        histDict["gginvmass"]["hist"].Fill(invMass)

      ###End of if TriggerPass
    ###End of event loop 

  f.Close()


foutName = outDir + outfilename + "_" + str(nBatchRun) + ".root"

print(foutName)

outFile = TFile.Open(foutName, "RECREATE")
outFile.cd()



###Write Hists
for hName, hDict in histDict.items():
  hDict["hist"].Write()


outFile.Close()
