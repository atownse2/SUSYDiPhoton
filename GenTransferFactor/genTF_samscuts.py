from utils import *

head = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/" 

#dType = "WGJets"
#dType = "DYJetsToLL_M-50"
dType = "TTJets"


tag = dataDict[dType]["tag"]
era = dataDict[dType]["era"]
HLT = dataDict[dType]["HLT"]


versiontag = "samscuts"

outDir = "hists/"
###I/O
ntuple_version =  "TreeMakerRandS_skimsv8"
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
  outDir = "/scratch365/atownse2/SUSYDiPhoton/GenTransferFactor/" + dType + "_" + era + "/"

if not any(filenames):
  print("No files in batch")
  quit()


print( "Running over " + str(len(filenames)) + " files")
###


###
outfilename = dType + "_" + era + "_" +  ntuple_version +  "_genTF_ns" 

if versiontag:
  outfilename += "_" + versiontag
if runbatch:
  outfilename += "_" + str(run)

outfilename += ".root"

print( "outfilename : " +  outfilename)
print( "outfiledir  : " +  outDir)
###


###


gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 




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
for filename in filenames:

  print("Opening File " + filename )
  f = TFile.Open(filename)
  t = f.Get('TreeMaker2/PreSelection')

  #print(t.ls())

  for ptRange in ptRanges:
    if ptRange in filename:
      ptTag = ptRange
      break
  
  if dType == "TTJets":
    ptTag = "all"

  tCounter = f.Get("tcounter")

  for i in range(tCounter.GetEntries()):
    hists["hEvents_"+ptTag].Fill(1)

  ##
  print("Starting Event Loop")
  for event in t:
    eventCount += 1

    if t.TriggerPass[HLT]==1:
      triggerPassCount += 1

      #if len(t.Photons) < 2: 
      #  continue

      EMpass = [ {"pt"  : t.Photons[i].Pt(),
                "4vec"  : t.Photons[i] , 
                "xseed" : bool(t.Photons_hasPixelSeed[i]),
                "barrel": t.Photons_isEB[i],
                "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 30 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

      #if len(EMpass) < 2:
      #  continue

      em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first

      #if not t.Photons_isEB[em[0]["index"]]  and not t.Photons_isEB[em[1]["index"]]:
      #  continue

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

          if deltaR(reco["4vec"], genEle) < 0.03 and pcdiff(reco["pt"], genEle.Pt()) < 20:
            
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
  print("Closing File")
  f.Close()
  

print("Event count = " + str(eventCount))


print("Writing hists to: "+outDir+outfilename)
outFile = TFile.Open(outDir+outfilename, "RECREATE")
outFile.cd()


###Write Hists
for h_name, hist in hists.items():
  hist.Write()

outFile.Close()
