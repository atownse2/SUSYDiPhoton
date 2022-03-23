from ROOT import TCanvas, TFile, TH1F, PyConfig, TLegend, TRatioPlot
from ROOT import gROOT

from collections import OrderedDict
from closureTestUtils import *

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

d = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/closure/"

dType = "WGJets"
#dType = "DYJetsToLL_M-50"
#dType = "TTJets_DiLept_TuneCUETP8M1"

#dType = "TTJets"

run = "Summer16v3"


fin = TFile.Open("hists/"+ dType + "_" + run + "_closure.root")



##Format hist Dict structure
hTypes = ["met", "met_lowBDT", "met_medBDT", "met_highBDT", "leadpt", "trailpt", "leadeta", "traileta"]

selections = ["scaled_eg", "gg"]


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
    h = TH1F(hName, hName, len(bins)-1, bins)

    for ptRange in ptRanges:
      hEvents = fin.Get("hEvents_"+ptRange)
      hName_ptRange = getHistName(sel, hType, ptRange)

      if sel == "scaled_eg":
        h_tmp_barrel = fin.Get(getHistName(sel, hType, ptRange, "barrel"))
        h_tmp_endcap = fin.Get(getHistName(sel, hType, ptRange, "endcap"))

        barrel_fakerate, barrel_fakerate_error = fakerateDict[run]["barrel"]
        endcap_fakerate, endcap_fakerate_error = fakerateDict[run]["endcap"]

        h_ptRange = TH1F(hName_ptRange, hName_ptRange, len(bins)-1, bins)
        h_ptRange.Add(h_tmp_barrel, h_tmp_endcap, barrel_fakerate, endcap_fakerate)
        
        for i in range(h_ptRange.GetNbinsX()):
          barrel_sys = barrel_fakerate_error*h_tmp_barrel.GetBinContent(i+1) 
          endcap_sys = endcap_fakerate_error*h_tmp_endcap.GetBinContent(i+1)
          stat = h_ptRange.GetBinError(i+1)
          
          bin_error = np.sqrt(barrel_sys**2 + endcap_sys**2 + stat**2)

          h_ptRange.SetBinError(i+1, bin_error)

      else:
        h_ptRange = fin.Get(hName_ptRange)

      h_ptRange.Scale(1/hEvents.GetEntries())
      h.Add(h_ptRange)

    histDict[hName] = h


##Drawing pt2 test
for hType in hTypes:
  c = TCanvas( hType, hType, 800, 800)

  c.SetTitle(hType)

  l = TLegend(0.6,0.8,0.9,0.9)

  h_gg_name = "gg_" + hType
  h_eg_name = "scaled_eg_" + hType

  h_gg = drawOverflow(histDict[h_gg_name])
  h_eg = drawOverflow(histDict[h_eg_name])

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

