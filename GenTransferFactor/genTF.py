import sys
from ROOT import gROOT, TFile
import numpy as np
from utils import dataDict, get_file_name, init_hists
from cuts import version_mapping


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


#version = "default"
cuts_version = "samscuts"


outDir = "hists/"
###I/O
#ntuple_version =  "TreeMakerRandS_skimsv8"
ntuple_version = "TreeMaker"

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
  outfilename = get_file_name(dType, era, ntuple_version, "genTF", cuts_version, run )
else:
  outfilename = get_file_name(dType, era, ntuple_version, "genTF", cuts_version)


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
    version_mapping[cuts_version](t, hists, HLT) 
  print("Closing File")
  f.Close()

print("Writing hists to: "+outDir+outfilename)
outFile = TFile.Open(outDir+outfilename, "RECREATE")
outFile.cd()


###Write Hists
for h_name, hist in hists.items():
  hist.Write()

outFile.Close()
