import sys
from ROOT import gROOT, TFile
import numpy as np
from utils import dataDict, get_file_name, get_hist_name, init_hists, deltaR, pcdiff


###
gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


head = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/" 

dType = "WGJets"
#dType = "DYJetsToLL_M-50"
#dType = "TTJets"


tag = dataDict[dType]["tag"]
era = dataDict[dType]["era"]
HLT = dataDict[dType]["HLT"]

version = "052522v4"

outDir = "hists/"
###I/O
ntuple_version =  "TreeMakerRandS_skimsv8"
#ntuple_version = "TreeMaker"

ntuple_location = "root://cmsxrootd.fnal.gov//store/group/lpcsusyphotons/" + ntuple_version + "/"

#txtfiledir = "/scratch365/atownse2/SUSYDiPhoton/condor/"
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
  outDir = "/scratch365/atownse2/SUSYDiPhoton/GenTransferFactor/" + dType + "_" + era + "/"

if not any(filenames):
  print("No files in batch")
  quit()


print( "Running over " + str(len(filenames)) + " files")
###


if runbatch:
  outfilename = get_file_name(dType, era, ntuple_version, "genTF", version, run )
else:
  outfilename = get_file_name(dType, era, ntuple_version, "genTF", version)


print( "outfilename : " +  outfilename)
print( "outfiledir  : " +  outDir)
###



### Initialize histograms
hists = init_hists(dType)

##Start Loop over files in batch
for filename in filenames:

  print("Opening File " + filename )
  f = TFile.Open(filename)
  t = f.Get('TreeMaker2/PreSelection')

  print("Starting Event Loop")
  for event in t:

    if t.TriggerPass[HLT] != 1:
      continue   

    EMpass = [ {"pt"  : t.Photons[i].Pt(),
                "eta" : t.Photons[i].Eta(),
                "4vec"  : t.Photons[i] , 
                "xseed" : bool(t.Photons_hasPixelSeed[i]),
                "barrel": t.Photons_isEB[i],
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 70 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

    if len(EMpass) < 2:
      continue

    em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first

    clean_jets = []

    for jet in [j for j in t.Jetsclean if j.Pt() > 30]:
      jetPass = True
      for e in em:
        if deltaR(jet, e["4vec"]) < 0.4:
          jetPass = False
          break
      if jetPass:
        clean_jets.append(jet)

    njets = len(clean_jets)

    vmult = t.NVtx

    #w = t.CrossSection*t.puWeight*1000*35.9
    w = 1


    for reco in em:
      if reco["barrel"]:
        region = "barrel"
      else:
        region = "endcap"

      genPar_name = "genJet"
      for iPar, genPar in enumerate(t.GenParticles):
        if genPar.Pt() < 10:
          continue

        if deltaR(reco["4vec"], genPar) > 0.2 or pcdiff(reco["pt"], genPar.Pt()) > 20:
          continue

        pdgid = abs(t.GenParticles_PdgId[iPar])

        if pdgid == 11:
          genPar_name = "genEle"
        elif pdgid == 13:
          genPar_name = "genMu"
        elif pdgid == 15:
          genPar_name = "genTau"
        elif pdgid == 22:
          genPar_name = "genPho"

        break

      if reco["xseed"]: #Electron
        #Fill Hists
        hists[get_hist_name(genPar_name, "recoEle", region, "pt" )].Fill(reco["4vec"].Pt() , w)
        hists[get_hist_name(genPar_name, "recoEle", region, "eta")].Fill(reco["4vec"].Eta(), w)
        hists[get_hist_name(genPar_name, "recoEle", region, "njets")].Fill(njets, w)
        hists[get_hist_name(genPar_name, "recoEle", region, "met")].Fill(t.MET, w)          
        hists[get_hist_name(genPar_name, "recoEle", region, "nEle")].Fill(0, w)

      else: #Reco Photon
        #Fill Hists
        hists[get_hist_name(genPar_name, "recoPho", region, "pt" )].Fill(reco["4vec"].Pt() , w)
        hists[get_hist_name(genPar_name, "recoPho", region, "eta")].Fill(reco["4vec"].Eta(), w)
        hists[get_hist_name(genPar_name, "recoPho", region, "njets")].Fill(njets, w)
        hists[get_hist_name(genPar_name, "recoPho", region, "met")].Fill(t.MET, w)
        hists[get_hist_name(genPar_name, "recoPho", region, "nPho")].Fill(0, w)

  print("Closing File")
  f.Close()

print("Writing hists to: "+outDir+outfilename)
outFile = TFile.Open(outDir+outfilename, "RECREATE")
outFile.cd()


###Write Hists
for h_name, hist in hists.items():
  hist.Write()

outFile.Close()
