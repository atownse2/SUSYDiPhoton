from ROOT  import RooRealVar, RooGenericPdf, RooArgList, RooVoigtian, RooAddPdf, RooDataHist, RooFit, TPaveText, TCanvas
from ROOT import gROOT

import sys
top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/"
sys.path.append(top_dir)

from analysis.logger import Logger

l = Logger()

gROOT.SetBatch() # don't pop up canvases
gROOT.SetStyle('Plain') 

def z_fit(dType, hist, plot=False):
  '''
  Fit ee/eg invariant mass hist with Voigt + Linear profile.
  Returns number of signal events.
  '''

  x = RooRealVar("x", "Invariant Mass", 65, 115)

  ###Model
  c0 = RooRealVar("c0", "Y intercept of linear bkg", 1, 0, 500)
  c1 = RooRealVar("c1", "Slope of linear bkg", 0.1, 0, 10)
  bkg = RooGenericPdf("bkg", "linear background", "c0 + c1*x", RooArgList(c0,c1,x))

  #Voigt Profile
  m = RooRealVar("m", "mean of Voigtian", 90, 86, 94)
  w = RooRealVar("w", "width of Voigtian", 0.1, 50)
  s = RooRealVar("s", "sigma of Voigtian", 0.1, 50)

  sig = RooVoigtian("voigt", "Voigtian", x, m, w, s)

  nsig = RooRealVar("nsig", "Number of signal events", 0, hist.GetEntries())
  nbkg = RooRealVar("nbkg", "Number of background events", 0, hist.GetEntries())

  model = RooAddPdf("model", "model", RooArgList(sig,bkg), RooArgList(nsig,nbkg))
  ###

  ###Fit
  hname = hist.GetName()
  frame = x.frame(RooFit.Title(hname))

  data = RooDataHist("dh","Invariant Mass", RooArgList(x), hist)
  fitResult = model.fitTo(data, RooFit.PrintLevel(4), RooFit.Save(True))

  if plot:
    data.plotOn(frame)  
    model.plotOn(frame)

    c = TCanvas(hname+"_c", dType+'_'+hname+'_c', 600, 600) 
    c.cd()

    frame.Draw()

    # Set the y-axis limits
    y_max = hist.GetMaximum() + hist.GetBinError(hist.GetMaximumBin())
    frame.SetMaximum(y_max*1.2)


    # Add a label to the plot
    label = TPaveText(0.7, 0.7, 0.9, 0.9, 'NDC')
    label.SetTextSize(0.03)  # Set the font size to 0.03
    label.SetTextAlign(33)   # Move the label to the top right

    label.AddText("N Signal = " + str(round(nsig.getVal(),2)))
    label.AddText("Err Sig  = " + str(round(nsig.getError(),2)))
    label.AddText("Fit Status: " + str(fitResult.status()))

    label.Draw()
    out = top_dir + "plots/fits/" + dType + '_' + hname + ".png"
    l.log(f"Saving plot to {out}", 2)
    c.SaveAs(out)

  return (nsig.getVal(), nsig.getError())