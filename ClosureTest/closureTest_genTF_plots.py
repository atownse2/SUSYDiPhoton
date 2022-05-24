from ROOT import TCanvas, TFile, TH1F, PyConfig, TLegend, TRatioPlot
from ROOT import gROOT

from collections import OrderedDict
from closureTestUtils import *

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 


fdir = "hists/"


#dType = "DYJetsToLL_M-50"
#run = "Summer16v3"

dType = "WGJets"
run = "Summer16v3"

#dType = "TTJets"
#run = "Summer16v3"



ntuple_version = "TreeMaker"
#ntuple_version = "TreeMakerRandS_skimsv8"

version_tag = "samscuts"
#version_tag = None


if version_tag is not None:
  fin_name  = fdir + dType + "_" + run + "_" + ntuple_version + "_closureTest_genTF_ns_" + version_tag +  ".root"
  fout_name = fdir + dType + "_" + run + "_" + ntuple_version + "_closureTest_genTF_"    + version_tag + ".root"
else:
  fin_name  = fdir + dType + "_" + run + "_" + ntuple_version + "_closureTest_genTF_ns.root"
  fout_name = fdir + dType + "_" + run + "_" + ntuple_version + "_closureTest_genTF.root"


fin = TFile.Open(fin_name)



genTransferFactor = {"barrel": (0.017383707625563526, 9.211294454590436e-05),
                     "endcap": (0.024732126508268247, 0.00017889773719544853)}



##Format hist Dict structure
hTypes = ["DiEMpt", "met", "met_lowBDT", "met_medBDT", "met_highBDT", "leadpt", "trailpt", "leadeta", "traileta"]

selections = ["eg", "gg"]

regions = ["barrel" , "endcap"]

tags = ["tagHpt" , "tagLpt"]

histDict = OrderedDict()





def getHistName(sel, hType, region = False, ptTag = False, bdt = False):
  hName = sel + "_" + hType + "_"
  if bdt:
    hName += bdt + "_"
  if region:
    hName += region + "_"
  if ptTag:
    hName += ptTag + "_"

  return hName

from closureTestUtils import varBins

for hType in hTypes:
  for varName, b in varBins.items():
    if varName in hType:
      bins = b

  for sel in selections:
    if sel == "eg":
      for region in regions:
        for tag in tags:
          hName = getHistName(sel, hType, region = region, ptTag = tag)
          histDict[hName] = fin.Get(hName)
    else:
      hName = getHistName(sel, hType)
      histDict[hName] = fin.Get(hName)










"""
#getHistName(sel, hType, ptRange, region = False, bdt = False):
# Load/Add/Scale Hists
for hType in hTypes:
  for varName, b in varBins.items():
    if varName in hType:
      bins = b
      break

  for sel in selections:
    hName = sel + "_" + hType

    if sel == "eg":
      for region in regions:
        for tag in tags:
          hName_region_tag = hName + "_" + region + "_" + tag
          h_region_tag = TH1F(hName_region_tag, hName_region_tag, len(bins)-1, bins)

          for ptRange in ptRanges:
            hEvents = fin.Get("hEvents_" + ptRange)
            hName_ptRange = getHistName("eg" , hType, ptRange, region, ptTag=tag)
            h_ptRange = fin.Get(hName_ptRange)

            h_ptRange.Scale(1/hEvents.GetEntries())
            h_region_tag.Add(h_ptRange)

          histDict[hName_region_tag] = h_region_tag

    else:
      h = TH1F(hName, hName, len(bins)-1, bins)

      for ptRange in ptRanges:
        hEvents = fin.Get("hEvents_" + ptRange)
        hName_ptRange = getHistName(sel, hType, ptRange)
        h_ptRange = fin.Get(hName_ptRange)
        h_ptRange.Scale(1/hEvents.GetEntries())
        h.Add(h_ptRange)

      histDict[hName] = h
"""


##Drawing pt2 test
for hType in hTypes:
  for varName, b in varBins.items():
    if varName in hType:
      bins = b
      break

  c = TCanvas( hType, hType, 800, 800)

  c.SetTitle(hType)

  l = TLegend(0.6,0.8,0.9,0.9)

  h_gg_name = getHistName("gg", hType)

  h_eg_name = getHistName("eg", hType)

  h_gg = drawOverflow(histDict[h_gg_name])

  h_eg = TH1F(h_eg_name, h_eg_name, len(bins)-1, bins)

  for region in regions:
    for tag in tags:
      tf = genTransferFactor[region][0]
      tf_plus_e = genTransferFactor[region][0] + genTransferFactor[region][1]

      hName_region_tag = getHistName("eg", hType, region, tag)
      h_region_tag = histDict[hName_region_tag]
      n_bins = h_region_tag.GetNbinsX()

      for ibin in range(n_bins+2): #Use Fakerate error to set error bars including under/overflow
        bin_tf = h_region_tag.GetBinContent(ibin)*tf
        bin_tf_plus_e = h_region_tag.GetBinContent(ibin)*tf_plus_e

        tf_sys_e = bin_tf_plus_e - bin_tf
        tf_stat = h_region_tag.GetBinError(ibin)*tf
        bin_error = np.sqrt(tf_sys_e**2 + tf_stat**2)

        h_region_tag.SetBinContent(ibin, bin_tf)
        h_region_tag.SetBinError(ibin, bin_error)

      h_eg.Add(h_region_tag)

  h_eg = drawOverflow(h_eg)

  h_gg.SetTitle(dType+hType)
  h_eg.SetTitle(dType+hType)

  #h_gg.SetStats(0)
  h_eg.SetStats(0)

  h_gg.SetLineColor(2)

  #maxEntries = max(h_gg.Integral(), h_eg.Integral())

  l.AddEntry(h_gg, "gg")
  l.AddEntry(h_eg, "Scaled eg")

  ggMax = h_gg.GetMaximum() + h_gg.GetBinError(h_gg.GetMaximumBin())
  egMax = h_eg.GetMaximum() + h_eg.GetBinError(h_eg.GetMaximumBin())

  yMax = max(ggMax,egMax)*1.0
  #yMax = 1

  hRatio = TRatioPlot(h_eg,h_gg) 

  hRatio.Draw()

  hRatio.GetUpperPad().cd()

  h_eg.Draw("hist E1 same")
  h_gg.Draw("hist E1 same")

  hRatio.GetUpperRefYaxis().SetRangeUser(0,yMax)
  hRatio.GetLowerRefGraph().SetMaximum(2)

  print(hType + " gg int: " + str(h_gg.Integral()))
  print(hType + " gg get Ent: " + str(h_gg.GetEntries()))

  l.Draw()

  print("Saving plots hopefully?")
  c.SaveAs("plots/"+dType+"_closureTest_genTF_" + version_tag + "_" + hType + ".png")
  #c.SaveAs(dType+"_closureTest_" + pType + ".png")



fout = TFile(fout_name, "RECREATE")

for hName, hist in histDict.items():
  hist.Write()

fin.Close()
fout.Close()

