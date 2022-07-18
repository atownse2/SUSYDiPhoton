from ROOT import TCanvas, TFile, TH1F, PyConfig, TLegend, TRatioPlot
from ROOT import gROOT

from utils import *
from hists import *

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

d = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src/closure/"

### Data Type

dType = "WGJets"
#dType = "DYJetsToLL_M-50"

#dType = "TTJets"

run = "Summer16v3"


###
version = "071822v2"

fin_dir = "hists/"
###I/O
#ntuple_version =  "TreeMakerRandS_skimsv8"
ntuple_version = "skimsforAustin"
fin_name = get_file_name(dType, run, ntuple_version, "cTest", version)


### Load and scale hists
load_hists(dType, fin_dir+fin_name)
scale_hists(dType)

##Drawing pt2 test
for hType in hTypes:
  bins = get_bins(hType)

  c = TCanvas( hType, hType, 800, 800)

  c.SetTitle(hType)

  l = TLegend(0.6,0.8,0.9,0.9)

  h_gg_name = "gg_" + hType
  h_eg_name = "eg_" + hType

  h_gg = scaled_hists[h_gg_name]

  h_eg = TH1F(h_eg_name, h_eg_name, len(bins)-1, bins)

  for region in regions:
    for tag in ptTags:
      tf = fakerateDict[run][region][tag][0]
      tf_plus_e = fakerateDict[run][region][tag][0] + fakerateDict[run][region][tag][1]

      hName_region_tag = h_eg_name + "_" + region + "_" + tag
      h_region_tag = scaled_hists[hName_region_tag]
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

  h_gg = drawOverflow(h_gg)
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
  c.SaveAs("plots/"+dType+"_"+run+"_"+ntuple_version+"_cTest_plots_"+version+"_"+hType+".png")
  #c.SaveAs(dType+"_closureTest_" + pType + ".png")



fout = TFile("hists/" + dType+ "_" + run +"_"+ntuple_version+"_cTest_plots_"+version+".root", "RECREATE")

for hName, hist in histDict.items():
  hist.Write()

fout.Close()

