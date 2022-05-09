from utils import * 

gStyle.SetOptStat(0)
gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

fdir = "hists/"


#dType = "DYJetsToLL_M-50"
#run = "Summer16v3"

#dType = "WGJets"
#run = "Summer16v3"

dType = "TTJets"
run = "Summer16v3"


version_tag = "samscuts"
#version_tag = None


if version_tag is not None:
  fin_name = fdir + dType + "_" + run + "_genTF_ns_" + version_tag +  ".root"
  fout_name = fdir + dType + "_" + run + "_genTF_" + version_tag + ".root"
else:
  fin_name = fdir + dType + "_" + run + "_genTF_ns.root"
  fout_name = fdir + dType + "_" + run + "_genTF.root"



print("File In: " + fin_name)
print("File Out: " + fout_name)


# Load hists
hists = load_scaled_hists(dType, fin_name)



### Write Hists
fout = TFile.Open(fout_name, "RECREATE")
fout.cd()

for hName, hist in hists.items():
  if hist is not None:
    hist.Write()

### Make plots
for genType in genTypes:
  for region in regions:
    for eVar in eVars:
      title = dType + "_" + genType + "_" + region + "_" + eVar
  
      if version_tag is not None:
        title +=  "_" + version_tag
 
      h_recoPho = hists[getHistName(genType, "recoPho" , region, eVar)]
      h_recoEle = hists[getHistName(genType, "recoEle" , region, eVar)]

      c = TCanvas(title+"_c", title+"_c", 1000, 800)
      c.SetLogy()
      l = TLegend(0.6,0.8,0.9,0.9)
      c.cd()

      h_recoPho.SetStats(0)
      h_recoEle.SetStats(0)

      h_recoEle.SetLineColor(2)

      l.AddEntry(h_recoPho, "recoPho")
      l.AddEntry(h_recoEle, "recoEle")

      recoPhoMax = h_recoPho.GetMaximum() + h_recoPho.GetBinError(h_recoPho.GetMaximumBin())
      recoEleMax = h_recoEle.GetMaximum() + h_recoEle.GetBinError(h_recoEle.GetMaximumBin())

      yMax = max(recoPhoMax,recoEleMax)*1.25

      hRatio = TRatioPlot(h_recoPho,h_recoEle) 
      hRatio.Draw()

      hRatio.GetUpperPad().cd()

      h_recoPho.Draw("hist E1 same")
      h_recoEle.Draw("hist E1 same")

      hRatio.GetUpperRefYaxis().SetRangeUser(1,yMax)
      #hRatio.GetLowerRefGraph().SetMaximum(0.05)

 

      c.Update()
      c.SaveAs("plots/" + title + ".png")


fout.Close()
