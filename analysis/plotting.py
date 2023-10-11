import sys

from ROOT import TCanvas, TRatioPlot, TLegend, TH1F, TH2F, gPad, gROOT

top_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/'
sys.path.append(top_dir)

plot_dir = top_dir + 'plots/'

def unfold_2d(h2d):
    x_bins = h2d.GetNbinsX()
    y_bins = h2d.GetNbinsY()
    hname = h2d.GetName() + "_unfolded"
    h1d = TH1F(hname, hname, x_bins*y_bins, 0, x_bins*y_bins)
    

def plot_egtf_closure(h_gg, h_gg_pred, 
                      gg_label=None, gg_pred_label=None, 
                      x_title=None, y_title=None, 
                      save=False):
    
    c = TCanvas( "c", "c", 800, 800)
    l = TLegend(0.6,0.8,0.9,0.9)

    # Unfold 2D hist into 1D



    h_eg.SetStats(0)
    h_gg.SetLineColor(2)
    h_eg.GetXaxis().SetTitle('hardMET')
    h_eg.SetTitle(title)

    l.AddEntry(h_gg, "gg")
    l.AddEntry(h_eg, "Scaled eg")

    hRatio = TRatioPlot(h_eg,h_gg) 
    hRatio.Draw()

    hRatio.GetUpperPad().cd()
    h_eg.Draw("hist E1 same")
    h_gg.Draw("hist E1 same")
    l.Draw()

    ggMax = h_gg.GetMaximum() + h_gg.GetBinError(h_gg.GetMaximumBin())
    egMax = h_eg.GetMaximum() + h_eg.GetBinError(h_eg.GetMaximumBin())
    yMax = max(ggMax,egMax)*1.2
    hRatio.GetUpperRefYaxis().SetRangeUser(0,yMax)
    #hRatio.GetLowerRefGraph().SetMaximum(2)

    c.SaveAs(plot_dir+title.replace(' ','_')+".png")

def plot_background_matrix(dType, hname, hist):
    c = TCanvas(hname + "_c", hname + "_c", 1000, 800)
    c.cd()

    x_title = "lead recoPhoton"
    y_title = "trail recoPhoton"

    hist.GetXaxis().LabelsOption("a")
    hist.GetYaxis().LabelsOption("a")

    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)

    hist.GetYaxis().SetTitleOffset(2.0)
    gPad.SetLeftMargin(0.15)

    hist.SetStats(False)
    hist.Draw("TEXT")

    fName = dType + "_"  + hname
    hist.SetTitle(fName)

    c.Update()
    c.SaveAs(plot_dir + fName+".png")


def compare_bdt_fakeBKG(title, h_EG, h_GJet):
    '''
    Compares bdt distribution of jet fakes to electron fakes.
    '''
    c = TCanvas(title, title, 1000, 800)

    h_EG.SetLineColor(1)
    h_GJet.SetLineColor(2)

    h_EG.Scale(1/h_EG.Integral())
    h_GJet.Scale(1/h_GJet.Integral())

    h_GJet.GetXaxis().SetTitle("BDT score")
    h_GJet.GetXaxis().SetRangeUser(-0.6,0.6)
    h_GJet.GetYaxis().SetTitle("Normalized Counts")
    h_GJet.GetYaxis().SetTitleOffset(1.25)
    h_GJet.SetTitle(title)

    hRatio = TRatioPlot(h_GJet,h_EG)
    hRatio.Draw()

    legend = TLegend(0.15,0.75,0.35,0.9)
    legend.AddEntry(h_EG, "genPho + genEle")
    legend.AddEntry(h_GJet, "genPho + genJet")
    legend.Draw()

    c.SaveAs(plot_dir+"{}.png".format(title.replace(' ', '_')))

def drawOverflow(h):
  nx = h.GetNbinsX() + 1
  xbins = array("d", [0]*(nx+1))
  for i in range(nx):
    xbins[i] = h.GetBinLowEdge(i+1)
  xbins[nx] = xbins[nx-1] + h.GetBinWidth(nx)
  tmpName = h.GetName() + "_Overflow"
  htmp = TH1F(tmpName, tmpName, nx, xbins)

  htmp.SetXTitle(h.GetXaxis().GetTitle())
  htmp.SetYTitle(h.GetYaxis().GetTitle())

  for j in range(nx):
    i = j+1
    htmp.Fill(htmp.GetBinCenter(i), h.GetBinContent(i))
    htmp.SetBinError(i, h.GetBinError(i))

  htmp.Fill(h.GetBinLowEdge(1)-1, h.GetBinContent(0))
  htmp.SetEntries(h.GetEntries())
  htmp.SetFillStyle(h.GetFillStyle())
  htmp.SetFillColor(h.GetFillColor())

  return htmp