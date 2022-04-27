from ROOT import TCanvas, TFile, TH1F, PyConfig, TLegend, TRatioPlot
from ROOT import gROOT

from collections import OrderedDict
from closureTestUtils import *

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

d = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/closure/"

dType = "WGJets"
#dType = "DYJetsToLL_M-50"

#dType = "TTJets"

run = "Summer16v3"


fin = TFile.Open("hists/"+ dType + "_" + run + "_closure.root")

if not fin:
  print("File not opened")
  quit()

##Format hist Dict structure
hTypes = ["DiEMpt", "met", "met_lowBDT", "met_medBDT", "met_highBDT", "leadpt", "trailpt", "leadeta", "traileta"]

selections = ["eg", "gg"]

regions = ["barrel" , "endcap"]

tags = ["tagHpt" , "tagLpt"]

histDict = OrderedDict()

ptRanges = ptRangeDict[dType]

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


##Drawing pt2 test
for hType in hTypes:
  for varName, b in varBins.items():
    if varName in hType:
      bins = b
      break

  c = TCanvas( hType, hType, 800, 800)

  c.SetTitle(hType)

  l = TLegend(0.6,0.8,0.9,0.9)

  h_gg_name = "gg_" + hType
  h_eg_name = "eg_" + hType

  h_gg = drawOverflow(histDict[h_gg_name])

  h_eg = TH1F(h_eg_name, h_eg_name, len(bins)-1, bins)

  for region in regions:
    for tag in tags:
      tf = 4*fakerateDict[run][region][tag][0]
      tf_plus_e = 4*fakerateDict[run][region][tag][0] + 4*fakerateDict[run][region][tag][1]

      hName_region_tag = h_eg_name + "_" + region + "_" + tag
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
  c.SaveAs("plots/"+dType+"_closureTest_" + hType + ".png")
  #c.SaveAs(dType+"_closureTest_" + pType + ".png")



fout = TFile("hists/" + dType+ "_" + run + "_closureTestPt2.root", "RECREATE")

for hName, hist in histDict.items():
  hist.Write()

fin.Close()
fout.Close()

