from ROOT import gROOT, TFile
import sys
import numpy as np

from utils import deltaR, dataDict, get_file_name
from hists import histDict, init_hists, fill_hists, write_hists, ptRangeDict


###
gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


head = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/" 

dType = "DYJetsToLL_M-50"

tag = dataDict[dType]["tag"]
era = dataDict[dType]["era"]
HLT = dataDict[dType]["HLT"]


script_name = "tf_MC"

version = "071822v1"

outDir = "hists/"
###I/O
#ntuple_version =  "TreeMakerRandS_skimsv8"
ntuple_version = "skimsforAustin"

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
  outDir = "/scratch365/atownse2/SUSYDiPhoton/TransferFactor/" + dType + "_" + era + "/"

if not any(filenames):
  print("No files in batch")
  quit()


print( "Running over " + str(len(filenames)) + " files")
###


if runbatch:
  outfilename = get_file_name(dType, era, ntuple_version, script_name, version, run )
else:
  outfilename = get_file_name(dType, era, ntuple_version, script_name, version)


print( "outfilename : " +  outfilename)
print( "outfiledir  : " +  outDir)
###

### Define hist dicts and bins
init_hists(dType, isMC=True)

ptRanges = ptRangeDict[dType]
###

##Start Loop over files in batch
for filename in filenames:
 
  f = TFile.Open(filename)
  print("Opening File " + filename )
  t = f.Get('TreeMaker2/PreSelection')

  #Parse Filename for "ptRange"
  for ptRange in ptRanges:
    if ptRange in filename:
      ptTag = ptRange
      break

  tCounter = f.Get("tcounter")

  for i in range(tCounter.GetEntries()):
    histDict["hEvents_"+ptTag].Fill(1)

  ##
  print("Starting Event Loop")
  for event in t:

    if t.TriggerPass[HLT]==1:

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

      #w = t.CrossSection*t.puWeight*1000*35.9

      lead  = em[0]
      trail = em[1]

      invmass = (lead["4vec"] + trail["4vec"]).M()

      if (lead["xseed"] and trail["xseed"]): #Passing probes ee

        if lead["barrel"] and trail["barrel"]:
          regions = ["BB"]
        elif lead["barrel"] != trail["barrel"]:
          regions = ["EB","BE"] #Not distinguishing between electrons
        else:
          regions = ["EE"]

        for region in regions:
          fill_hists( "ee" , region, invmass, ptRange=ptTag) ##TODO add ability to look at fakerate as fn of other kinematic variables


      elif (lead["xseed"] != trail["xseed"]): #Failing probes eg
        electron = lead if lead["xseed"] else trail
        photon = trail if lead["xseed"] else lead

        region = ("B" if electron["barrel"] else "E") + ("B" if photon["barrel"] else "E")

        fill_hists("eg" , region, invmass, ptRange=ptTag)

      ###End of if TriggerPass
    ###End of event loop 
  f.Close()

write_hists(outDir+outfilename)
