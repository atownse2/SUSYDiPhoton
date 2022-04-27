from ROOT import TCanvas, TFile, TH1F, TH1I, PyConfig
from ROOT import gROOT
import sys
import numpy as np
from utils import dR

dType = "DYJetsToLL_M-50"
run = "Summer16v3"
rtag = "HT"
HLT = 21

#txtfiledir = "/scratch365/atownse2/SUSYDiPhoton/condor/"
txtfiledir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/samples/"

#outDir = "/scratch365/atownse2/SUSYDiPhoton/fakerate/" + dType + "_" + run + "/"
outDir = "hists/"

outfilename = dType + "_" + run + "_fakerate_ns"
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
from init_hists_mc import *
init_hists(dType)


genEle_recoPho_barrel = TH1I("genEle_recoPho_barrel" , "genEle_recoPho_barrel" , 1 , 0 , 1)
genEle_recoEle_barrel = TH1I("genEle_recoEle_barrel" , "genEle_recoEle_barrel" , 1 , 0 , 1)

genEle_barrel = TH1I("genEle_barrel" , "genEle_barrel" , 1 , 0 , 1)

genEle_recoPho_endcap = TH1I("genEle_recoPho_endcap" , "genEle_recoPho_endcap" , 1 , 0 , 1)
genEle_recoEle_endcap = TH1I("genEle_recoEle_endcap" , "genEle_recoEle_endcap" , 1 , 0 , 1)

genEle_endcap = TH1I("genEle_endcap" , "genEle_endcap" , 1 , 0 , 1)

genJet_recoPho_barrel = TH1I("genJet_recoPho_barrel" , "genJet_recoPho_barrel" , 1 , 0 , 1)
genJet_recoEle_barrel = TH1I("genJet_recoEle_barrel" , "genJet_recoEle_barrel" , 1 , 0 , 1)

genJet_barrel = TH1I("genJet_barrel" , "genJet_barrel" , 1 , 0 , 1)

genJet_recoPho_endcap = TH1I("genJet_recoPho_endcap" , "genJet_recoPho_endcap" , 1 , 0 , 1)
genJet_recoEle_endcap = TH1I("genJet_recoEle_endcap" , "genJet_recoEle_endcap" , 1 , 0 , 1)

genJet_endcap = TH1I("genJet_endcap" , "genJet_endcap" , 1 , 0 , 1)

genFakerateHists = [genEle_recoPho_barrel, genEle_recoEle_barrel, genEle_barrel, genEle_recoPho_endcap, genEle_recoEle_endcap, genEle_endcap, genJet_recoPho_barrel, genJet_recoEle_barrel, genJet_barrel, genJet_recoPho_endcap, genJet_recoEle_endcap, genJet_endcap]


ptRanges = ptRangeDict[dType]
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
    histDict["hEvents_"+ptTag]["hist"].Fill(1)

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

      nJets = len(jets)

      vmult = t.NVtx

      MET = t.MET

      w = t.CrossSection*t.puWeight*1000*35.9

      ## Generator level fakerate
      for genE in t.GenElectrons:
        if genE.Pt() < 70:
          continue 

        if abs(genE.Eta()) < 1.48:
          genEle_barrel.Fill(0)
          for recoE in em:
            if dR(recoE["4vec"], genE) < 0.4 and recoE["barrel"]:
              if recoE["xseed"]:
                genEle_recoEle_barrel.Fill(0)
              else:
                genEle_recoPho_barrel.Fill(0)
              break          

        else:
          genEle_endcap.Fill(0)
          for recoE in em:
            if dR(recoE["4vec"], genE) < 0.4 and not recoE["barrel"]:
              if recoE["xseed"]:
                genEle_recoEle_endcap.Fill(0)
              else:
                genEle_recoPho_endcap.Fill(0)              
              break 

      for genJ in t.GenJets:
        if genJ.Pt() < 70:
          continue

        ###Reject actual electrons and photons
        egMatch = False
        for iGen, genP in enumerate(t.GenParticles):
          if genP.Pt() < 70:
            continue
          if (t.GenParticles_PdgId[iGen] != 22 and abs(t.GenParticles_PdgId[iGen]) != 11):
            continue

          if dR(genJ, genP) < 0.4:
            egMatch = True
            break

        if egMatch: 
          continue

        ###Hopefully only hadronic jets now

        if abs(genJ.Eta()) < 1.48:
          genJet_barrel.Fill(0)
          for recoE in em:
            if dR(recoE["4vec"], genP) < 0.4 and recoE["barrel"]:
              if recoE["xseed"]:
                genJet_recoEle_barrel.Fill(0)
              else:
                genJet_recoPho_barrel.Fill(0)
              break          

        else:
          genJet_endcap.Fill(0)
          for recoE in em:
            if dR(recoE["4vec"], genP) < 0.4 and not recoE["barrel"]:
              if recoE["xseed"]:
                genJet_recoEle_endcap.Fill(0)
              else:
                genJet_recoPho_endcap.Fill(0)              
              break 

      ##

      lead  = em[0]
      trail = em[1]

      invMass = (lead["4vec"] + trail["4vec"]).M()

      if (lead["xseed"] and trail["xseed"]): #Passing probes

        fill_hist(ptTag, "ee", "pt" , trail["pt"]    , invMass, w, tag = "tagHpt") #Tag is H pt (lead), probe is L pt (trail)
        fill_hist(ptTag, "ee", "eta", trail["barrel"], invMass, w, tag = "tagHpt")
  
        fill_hist(ptTag, "ee", "pt" , lead["pt"]    , invMass, w, tag = "tagLpt") #Tag is L pt (trail), probe is H pt (lead)      
        fill_hist(ptTag, "ee", "eta", lead["barrel"], invMass, w, tag = "tagLpt")

        fill_hist(ptTag, "ee" , "jnum" , nJets, invMass, w)
        fill_hist(ptTag, "ee" , "vmult", vmult, invMass, w)
        fill_hist(ptTag, "ee" , "met"  , MET  , invMass, w)

      elif (lead["xseed"] != trail["xseed"]): #Failing probes

        if lead["xseed"]:
          tag = lead
          probe = trail
        else:
          tag = trail
          probe = lead

        if tag["pt"] > probe["pt"]: #tagHpt
          fill_hist(ptTag, "eg", "pt" , probe["pt"]    , invMass, w, tag = "tagHpt")
          fill_hist(ptTag, "eg", "eta", probe["barrel"], invMass, w, tag = "tagHpt")

        else: #tagLpt
          fill_hist(ptTag, "eg", "pt" , probe["pt"]    , invMass, w, tag = "tagLpt")
          fill_hist(ptTag, "eg", "eta", probe["barrel"], invMass, w, tag = "tagLpt")

        fill_hist(ptTag, "eg", "jnum" , nJets, invMass, w)
        fill_hist(ptTag, "eg", "vmult", vmult, invMass, w)
        fill_hist(ptTag, "eg", "met"  , MET  , invMass, w)

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
for hName, hDict in histDict.items():
  hDict["hist"].Write()

for hist in genFakerateHists:
  hist.Write() 

outFile.Close()

"""
      tags = [i for i in em if i["xseed"]]

      for tag in tags:
        for probe in em:
          if dR(tag["4vec"], probe["4vec"]) < 0.4:
            continue

          invMass = (tag["4vec"] + probe["4vec"]).M()

          if probe["xseed"]: #Passing probes
          
            if tag["pt"] > probe["pt"]: #tagHpt
              fill_hist(ptTag, "ee", "pt" , probe["pt"]    , invMass, w, tag = "tagHpt")
              fill_hist(ptTag, "ee", "eta", probe["barrel"], invMass, w, tag = "tagHpt")

            else: #tagLpt            
              fill_hist(ptTag, "ee", "pt" , probe["pt"]    , invMass, w, tag = "tagLpt")       
              fill_hist(ptTag, "ee", "eta", probe["barrel"], invMass, w, tag = "tagLpt")

            fill_hist(ptTag, "ee", "pt" , probe["pt"]    , invMass, w, tag = "all") 
            fill_hist(ptTag, "ee", "eta", probe["barrel"], invMass, w, tag = "all") #Using this to calculate transfer factor

            fill_hist(ptTag, "ee" , "jnum" , nJets, invMass, w)
            fill_hist(ptTag, "ee" , "vmult", vmult, invMass, w)
            fill_hist(ptTag, "ee" , "met"  , MET  , invMass, w)

          else: #Failing probes
            
            if tag["pt"] > probe["pt"]: #tagHpt
              fill_hist(ptTag, "eg", "pt" , probe["pt"]    , invMass, w, tag = "tagHpt")
              fill_hist(ptTag, "eg", "eta", probe["barrel"], invMass, w, tag = "tagHpt")

            else: #tagLpt
              fill_hist(ptTag, "eg", "pt" , probe["pt"]    , invMass, w, tag = "tagHpt")
              fill_hist(ptTag, "eg", "eta", probe["barrel"], invMass, w, tag = "tagHpt")

            fill_hist(ptTag, "eg", "pt" , probe["pt"]    , invMass, w, tag = "all") 
            fill_hist(ptTag, "eg", "eta", probe["barrel"], invMass, w, tag = "all") #Using this to calculate transfer factor 

            fill_hist(ptTag, "eg", "jnum" , nJets, invMass, w)
            fill_hist(ptTag, "eg", "vmult", vmult, invMass, w)
            fill_hist(ptTag, "eg", "met"  , MET  , invMass, w)
"""
