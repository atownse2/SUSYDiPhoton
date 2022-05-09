import numpy as np
import pickle
import sys

from ROOT import TH1F, TFile, Math, TCanvas, TLegend, TRatioPlot, gStyle, gROOT 
from collections import OrderedDict
from array import array

dataDict = {"WGJets": {"era": "Summer16v3",
                       "tag": "PtG",
                       "HLT": 21},
            "TTJets": {"era": "Summer16v3",
                       "tag": "",
                       "HLT": 21},
            "DYJetsToLL_M-50":{"era": "Summer16v3",
                       "tag": "HT",
                       "HLT": 21}}


### Define hist dicts and bins
ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets": ["all"]}

varBins = { "pt"    : np.linspace(0,400, num = 5),
            "eta"   : np.linspace(-2.4, 2.4, num=24),
            "njets"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250]}

genTypes = ["genEle"]

recoTypes = ["recoPho", "recoEle"]

regions = ["barrel" , "endcap"]

eVars = ["pt", "eta" , "njets" , "met"]

def getHistName(genType, recoType, region, var, ptRange = False):
  histName = genType + "_" + recoType + "_" + region + "_" + var

  if ptRange:
    histName += "_" + ptRange

  return histName


def init_hists(dType):
  ptRanges = ptRangeDict[dType]

  hists = OrderedDict()

  for ptRange in ptRanges:
    hists["hEvents_"+ptRange] = TH1F("hEvents_"+ptRange, "hEvents_"+ptRange, 1, 0, 1 ) 

    for genType in genTypes:
      for recoType in recoTypes:
        for region in regions:
          for eVar in eVars:
            h_name = getHistName(ptRange, genType, recoType, region, eVar)
            h_bins = array("d", varBins[eVar])
            hists[h_name] = TH1F( h_name, h_name, len(h_bins)-1, h_bins)

  return hists

def load_scaled_hists(dType, filename):

  f = TFile.Open(filename)

  ptRanges = ptRangeDict[dType]

  hists = OrderedDict()

  hEvents = { ptRange : f.Get("hEvents_"+ptRange).GetEntries() for ptRange in ptRanges }

  for genType in genTypes:
    for recoType in recoTypes:
      for region in regions:
        for eVar in eVars:
          h_name = getHistName(genType, recoType, region, eVar)
          h_bins = array("d", varBins[eVar])
          h = TH1F( h_name, h_name, len(h_bins)-1, h_bins)
          h.SetDirectory(0) #Write to RAM
          for ptRange in ptRanges:
            h_name_pt = getHistName(genType, recoType, region, eVar, ptRange)
            h_pt = f.Get(h_name_pt)
            print( h_name_pt )
            h_pt.Scale(2.5/hEvents[ptRange])
            h.Add(h_pt) 
          hists[h_name] = h
  f.Close()
  return hists


def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)




def transferFactor(eeSig, eeSigError, egSig, egSigError):
  fakerate = egSig/eeSig

  EPeg = 1/(eeSig)
  EPee = egSig/((eeSig**2))

  error = np.sqrt((EPeg*egSigError)**2 + (EPee*eeSigError)**2)

  return fakerate, error

def transferFactor2(eeSig, eeSigError, egSig, egSigError):
  fakerate = 0.5*egSig/eeSig

  EPeg = 1/(2*eeSig)
  EPee = egSig/(2*(eeSig**2))

  error = np.sqrt((EPeg*egSigError)**2 + (EPee*eeSigError)**2)

  return fakerate, error



### Utility functions
def deltaR(photon,jet):

  return Math.VectorUtil.DeltaR(photon, jet)
  """
  deltaPhi = abs((photon.Phi() - jet.Phi()))
  if deltaPhi > np.pi : 
    deltaPhi = 2*np.pi - deltaPhi

  deltaPhi2 = deltaPhi**2
  deltaEta2 = (photon.Eta() - jet.Eta())**2
  deltaR = np.sqrt(deltaPhi2 + deltaEta2)
  return deltaR
  """

def getBDTbin(bdt_val):
  if bdt_val < -0.13:
    return "lowBDT"
  elif bdt_val < 0.03:
    return "medBDT"
  else:
    return "highBDT"

### Fitting Invariant mass plots
from ROOT import *

def fitHist(hmass):
   
    hName = hmass.GetName() 

    print("Fitting " + hName)

    hmass.SetLineWidth(2)
    hmass.SetLineColor(kBlack)
    hmass.SetMarkerColor(kBlack)
    hmass.GetXaxis().SetTitle('m(e,#gamma) [GeV]')
    if not hmass.GetSumw2N(): hmass.Sumw2()    
    xax = hmass.GetXaxis()

    lowedge, highedge = 65,115

    x = RooRealVar("x", "Invariant Mass", lowedge, highedge)

    data = RooDataHist("dh","Invariant Mass", RooArgList(x), hmass)

    frame = x.frame(RooFit.Title(hName))
    data.plotOn(frame)  

    m = RooRealVar("m", "mean of Voigtian", 90, 88, 92)
    w = RooRealVar("w", "width of Voigtian", 0, 50)
    s = RooRealVar("s", "sigma of Voigtian", 0, 50)

    c0 = RooRealVar("c0", "Y intercept of linear bkg", 1, 0, 500)
    c1 = RooRealVar("c1", "Slope of linear bkg", 0.1, 0, 10)
  
    bkg = RooGenericPdf("bkg", "linear background", "c0 + c1*x", RooArgList(c0,c1,x))

    #fsig = RooRealVar("fsig", "fraction of signal", 0.5, 0, 1)

    sig = RooVoigtian("voigt", "Voigtian", x, m, w, s)

    nsig = RooRealVar("nsig", "Number of signal events", 0, hmass.GetEntries())
    nbkg = RooRealVar("nbkg", "Number of background events", 0, hmass.GetEntries())

    model = RooAddPdf("model", "model", RooArgList(sig,bkg), RooArgList(nsig,nbkg))

    fitResult = model.fitTo(data, RooFit.PrintLevel(-1), RooFit.Save(True))
    model.plotOn(frame)

    print("Number of Expected Signal Events = " + str(nsig.getVal()))
    print("Error  on Expected Signal Events = " + str(nsig.getError()))

    print("Number of Expected Background Events = " + str(nbkg.getVal()))

    print(" ")

    return nsig.getVal() , nsig.getError(), fitResult.status(), frame
