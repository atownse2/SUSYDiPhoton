import sys

from ROOT import gROOT, TFile
import numpy as np
from utils import dataDict, get_file_name, init_hists, genPar_map, passID, deltaR, pcdiff


###
gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


head = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/" 

#dType = "WGJets"
#dType = "DYJetsToLL_M-50"
dType = "TTJets"


tag = dataDict[dType]["tag"]
era = dataDict[dType]["era"]
HLT = dataDict[dType]["HLT"]

version = "060822v2"

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
  outDir = "/scratch365/atownse2/SUSYDiPhoton/FakeBackground/" + dType + "_" + era + "/"

if not any(filenames):
  print("No files in batch")
  quit()


print( "Running over " + str(len(filenames)) + " files")
###


if runbatch:
  outfilename = get_file_name(dType, era, ntuple_version, "fakeBKG", version, run )
else:
  outfilename = get_file_name(dType, era, ntuple_version, "fakeBKG", version)


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

    if t.IsRandS != 0:
      continue

    if t.TriggerPass[HLT] != 1:
      continue   

    EMpass = [ {"pt"  : t.Photons[i].Pt(),
                "eta" : t.Photons[i].Eta(),
                "4vec"  : t.Photons[i] , 
                "xseed" : bool(t.Photons_hasPixelSeed[i]),
                "barrel": t.Photons_isEB[i],
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 70 and abs(t.Photons[i].Eta()) < 2.4 and bool(t.Photons_fullID[i]) ] #passID(era, "tight", t, i)


    if len(EMpass) < 2:
      continue

    em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first
    
    lead = em[0]
    trail = em[1]

    if not any([lead["barrel"], trail["barrel"]]):
      continue

    #Count events
    hists["num_events"].Fill(0)
    
    genParticles = []
    for iPar, genPar in enumerate(t.GenParticles):
      if genPar.Pt() < 10:
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
      else:
        continue

      genParticles.append({ "genPar_name": genPar_name,
                            "genPar"     : genPar})


    for reco in [lead, trail]: #Match to genParticles
      
      for gen in genParticles:
        gen["dR"] = deltaR(reco["4vec"], gen["genPar"])
        gen["ptDiff"] = pcdiff(reco["pt"], gen["genPar"].Pt())

      passdR = [ g for g in genParticles if g["dR"] < 0.4 ]
 
      reco["gen_passdR"] = passdR

      if len(passdR) == 1:
        reco["gen_match"] = passdR[0]
        reco["gen_match_name"] = passdR[0]["genPar_name"]

      elif len(passdR) > 1:
        gen_pt_sort = sorted( passdR, key = lambda i: i["ptDiff"] ) ##Associate match with closest in pt
        reco["gen_match"] = gen_pt_sort[0]
        reco["gen_match_name"] = gen_pt_sort[0]["genPar_name"]
      else:
        reco["gen_match_name"] = "genJet"
        
        #For no matches look at dR and hadTowOverEM dist recoType + "_" + region + "_noMatch_" + genType + "_" + eVar
        region = "barrel" if reco["barrel"] else "endcap"
        recoType = "recoEle" if reco["xseed"] else "recoPho"
        gen_dR_sort = sorted(genParticles, key = lambda i: i["dR"])

        hists[recoType + "_" + region + "_noMatch_" + gen_dR_sort[0]["genPar_name"] + "_dR"].Fill(gen_dR_sort[0]["dR"])
        hists[recoType + "_" + region + "_noMatch_" + gen_dR_sort[0]["genPar_name"] + "_hadTowOverEM"].Fill(t.Photons_hadTowOverEM[reco["index"]])



    g_lead  = lead["gen_match_name"] 
    g_trail = trail["gen_match_name"] 

    if lead["barrel"] and trail["barrel"]:
      region2 = "bb"
    elif lead["barrel"] and not trail["barrel"]:
      region2 = "be"
    elif not lead["barrel"] and trail["barrel"]:
      region2 = "eb"
    else:
      print("Both objects are in the endcap")
      quit()

    if lead["xseed"] and not trail["xseed"]:
      fstate = "eg"
    elif not lead["xseed"] and trail["xseed"]:
      fstate = "ge"
    elif not lead["xseed"] and not trail["xseed"]:
      fstate = "gg"
    else:
      fstate = "ee"

    hists[fstate + "_" + region2 + "_matrix"].Fill(g_lead, g_trail, 1)


  print("Closing File")
  f.Close()

print("Writing hists to: "+outDir+outfilename)
outFile = TFile.Open(outDir+outfilename, "RECREATE")
outFile.cd()


###Write Hists
for h_name, hist in hists.items():
  hist.Write()

outFile.Close()
