from ROOT import *
from utils import *
import os, sys, math
import numpy as np
import pickle
gStyle.SetOptStat(0)
gROOT.SetBatch(1)

fDir = "hists/"

dType = "DoubleEG"
run = "2016"

#dType = "DoubleEG"
#run  = "2017"

#dType = "EGamma"
#run = "2018"

dName = dType + "_Run" + run

fname = dName + "_fakerate.root"

infile = TFile(fDir+fname)

#infile.ls()
keys = infile.GetListOfKeys()

dataDict = { "eeinvmass":{"hist":infile.Get("eeinvmass")},
             "eginvmass":{"hist":infile.Get("eginvmass")},
             "gginvmass":{"hist":infile.Get("gginvmass")}}


def fitHist(hmass):

    hmass.Rebin()  
 
    hName = hmass.GetName() 


    print("Fitting " + hName)

    hmass.SetLineWidth(2)
    hmass.SetLineColor(kBlack)
    hmass.SetMarkerColor(kBlack)
    hmass.GetXaxis().SetTitle('m(e,#gamma) [GeV]')
    if not hmass.GetSumw2N(): hmass.Sumw2()    
    xax = hmass.GetXaxis()

    lowedge, highedge = 70,110
    
    leg = mklegend(x1=.52, y1=.55, x2=.88, y2=.78, color=kWhite)

    x = RooRealVar("x", "Invariant Mass", lowedge, highedge)

    data = RooDataHist("dh","Invariant Mass", RooArgList(x), hmass)

    frame = x.frame(RooFit.Title(hName))
    data.plotOn(frame)  

    m = RooRealVar("m", "mean of Voigtian", 90, 88, 92)
    w = RooRealVar("w", "width of Voigtian", 2.7, 0, 8)
    s = RooRealVar("s", "sigma of Voigtian", 1.7, 0, 5)

    c0 = RooRealVar("c0", "Y intercept of linear bkg", 1, 0, 500)
    c1 = RooRealVar("c1", "Slope of linear bkg", 0.1, 0, 10)
  
    bkg = RooGenericPdf("bkg", "linear background", "c0 + c1*x", RooArgList(c0,c1,x))

    sig = RooVoigtian("voigt", "Voigtian", x, m, w, s)

    nsig = RooRealVar("nsig", "Number of signal events", 0, hmass.GetEntries())
    nbkg = RooRealVar("nbkg", "Number of background events", 0, hmass.GetEntries())

    model = RooAddPdf("model", "model", RooArgList(sig,bkg), RooArgList(nsig,nbkg))

    fitResult = model.fitTo(data, RooFit.PrintLevel(-1), RooFit.Save(True))
    model.plotOn(frame)

    print("Number of Expected Signal Events = " + str(nsig.getVal()))
    print("Error  on Expected Signal Events = " + str(nsig.getError()))

    print("Number of Expected Background Events = " + str(nbkg.getVal()))

    print("Mean Voigt " + str(m.getVal()))
    print("Width Voigt " + str(w.getVal()))
    print("Sigma Voigt " + str(s.getVal()))

    print(" ")

    cName = dType + "_Run" + run + "_" + hName
    c = TCanvas("c", cName, 600, 600)
  
    c.cd()

    frame.Draw()

    label = TPaveText( 0.9, 0.7, 0.9, 0.80, 'NDC') 
    label.AddText("N Signal = " + str(nsig.getVal()))
    label.AddText("Err Sig  = " + str(nsig.getError()))
    label.Draw()

    c.SaveAs("fits/" + cName + ".png")
    c.Close()

    
    return data, model, sig, bkg, nsig, nbkg


### Fit ee,eg,gginvmass hists
data_ee, model_ee, sig_ee, bkg_ee, nsig_ee, nbkg_ee = fitHist(dataDict["eeinvmass"]["hist"])
data_eg, model_eg, sig_eg, bkg_eg, nsig_eg, nbkg_eg = fitHist(dataDict["eginvmass"]["hist"])
data_gg, model_gg, sig_gg, bkg_gg, nsig_gg, nbkg_gg = fitHist(dataDict["gginvmass"]["hist"])


### Calculate the fakerate
fakerate = nsig_eg.getVal()/nsig_ee.getVal()
print("Fakerate for " + dType + " run " + run + " is : " + str(np.round(fakerate, decimals = 3) ))

### Plot eg peak fit on gg data

model = RooAddPdf("model", "model", RooArgList(sig_eg,bkg_gg), RooArgList(nsig_eg,nbkg_gg))

c = TCanvas("c" , "Scaled eg peak on gg data", 600, 600)
c.cd()

lowedge, highedge = 70,110
x = RooRealVar("x", "Invariant Mass", lowedge, highedge)
frame = x.frame(RooFit.Title("Invariant Mass"))

#data_eg.plotOn(frame)
#model.plotOn(frame)

#c.SaveAs("test.png")
