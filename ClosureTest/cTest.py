from ROOT import TCanvas, TFile, TH1F, PyConfig
from ROOT import gROOT

from utils import *
from hists import *


head = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/" 



gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


dType = "WGJets"
era = "Summer16v3"
tag = "PtG"
HLT = 21

#dType = "TTJets"
#run = "Summer16v3"
#ptTag = "genMET-150"
#HLT = 21





version = "071822v2"

outdir = "hists/"
###I/O
#ntuple_version =  "TreeMakerRandS_skimsv8"
ntuple_version = "skimsforAustin"

redir = "ndcms.crc.nd.edu"
#redir = "cmsxrootd.fnal.gov"
ntuple_location = "root://"+redir+"//store/group/lpcsusyphotons/" + ntuple_version + "/"

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
  outdir = "/scratch365/atownse2/SUSYDiPhoton/ClosureTest/" + dType + "_" + era + "/"

if not any(filenames):
  print("No files in batch")
  quit()


print( "Running over " + str(len(filenames)) + " files")
###


if runbatch:
  outfilename = get_file_name(dType, era, ntuple_version, "cTest", version, run )
else:
  outfilename = get_file_name(dType, era, ntuple_version, "cTest", version)


print( "outfilename : " +  outfilename)
print( "outfiledir  : " +  outdir)
###


### Initialize histograms
hists = init_hists(dType)
ptRanges = ptRangeDict[dType]

##Start Loop over files in batch
for filename in filenames:
 
  f = TFile.Open(filename)
  print("Opening File " + filename )
  t = f.Get('TreeMaker2/PreSelection')

  for pt in ptRanges:
    if pt in filename:
      ptRange = pt
      break

  if not ptRange:
    print("Tag not found:" + filename)
    quit()

  tCounter = f.Get("tcounter")

  for i in range(tCounter.GetEntries()):
    histDict["hEvents_"+ptRange].Fill(1)

  ##
  print("Starting Event Loop")
  for event in t:

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

      ### Require to be matched to either photon or electron
      lead["fail"] = True
      trail["fail"] = True

      for iPar, genPar in enumerate(t.GenParticles):
        if genPar.Pt() < 10:
          continue

        pdgid = abs(t.GenParticles_PdgId[iPar])

        if pdgid != 11 and pdgid != 22:
          continue

        for reco in [lead, trail]:
          if deltaR(reco["4vec"], genPar) < 0.4:
            reco["fail"] = False

      if (lead["fail"] or trail["fail"]):
        continue

      ###


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

        histDict[getHistName("eg", "met", region = region, ptTag = tag, ptRange = ptRange)].Fill(MET, w)
        histDict[getHistName("eg", "met", region = region, ptTag = tag, ptRange = ptRange, bdt = BDT)].Fill(MET, w)

        histDict[getHistName("eg", "BDT", region = region, ptTag = tag, ptRange = ptRange)].Fill(BDT_val, w)

        histDict[getHistName("eg", "pt_lead", region = region, ptTag = tag, ptRange = ptRange)].Fill(lead["pt"], w) 
        histDict[getHistName("eg", "pt_trail", region = region, ptTag = tag, ptRange = ptRange)].Fill(trail["pt"], w)
        histDict[getHistName("eg", "eta_lead", region = region, ptTag = tag, ptRange = ptRange)].Fill(abs(lead["eta"]), w)
        histDict[getHistName("eg", "eta_trail", region = region, ptTag = tag, ptRange = ptRange)].Fill(abs(trail["eta"]), w)

        DiEMpt = (lead["4vec"] + trail["4vec"]).M()

        histDict[getHistName("eg", "DiEMpt", region = region, ptTag = tag, ptRange = ptRange)].Fill(DiEMpt, w)


      elif not lead["xseed"] and not trail["xseed"]: # gamma + gamma

        histDict[getHistName("gg", "met", ptRange = ptRange)].Fill(MET, w)
        histDict[getHistName("gg", "met", ptRange = ptRange, bdt = BDT)].Fill(MET, w)

        histDict[getHistName("gg", "BDT", ptRange = ptRange)].Fill(BDT_val, w)

        histDict[getHistName("gg", "pt_lead", ptRange = ptRange)].Fill(lead["pt"], w)
        histDict[getHistName("gg", "pt_trail", ptRange = ptRange)].Fill(trail["pt"], w)
        histDict[getHistName("gg", "eta_lead", ptRange = ptRange)].Fill(abs(lead["eta"]), w)
        histDict[getHistName("gg", "eta_trail", ptRange = ptRange)].Fill(abs(trail["eta"]), w)

        DiEMpt = (lead["4vec"] + trail["4vec"]).M()

        histDict[getHistName("gg", "DiEMpt", ptRange = ptRange)].Fill(DiEMpt, w)


      ###End of if TriggerPass
    ###End of event loop 

  f.Close()


###Write Hists
write_hists(outdir + outfilename)

