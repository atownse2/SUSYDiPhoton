import numpy as np
import pickle


from ROOT import TH1F, TH2I, TFile, Math, TCanvas, TLegend, TRatioPlot, gStyle, gROOT 
import ROOT as ROOT
from collections import OrderedDict
from array import array



dataDict = {"WGJets": {"era": "Summer16v3",
                       "tag": "",
                       "HLT": 21},
            "TTJets": {"era": "Summer16v3",
                       "tag": "",
                       "HLT": 21},
            "DYJetsToLL_M-50":{"era": "Summer16v3",
                       "tag": "HT",
                       "HLT": 21}}

def get_file_name(dType, era, ntuple_version, script_version, git_version, nbatch=None):
  f = dType + "_" + era + "_" + ntuple_version + "_" + script_version + "_" + git_version

  if nbatch is not None:
    f += "_" + str(nbatch)

  return f + ".root" 

### Define hist dicts and bins

genPar_map = { "genPho" : 0,
               "genEle" : 1,
               "genMu"  : 2,
               "genTau" : 3,
               "genJet" : 4}

varBins = { "pt"    : np.linspace(0,400, num = 32),
            "eta"   : np.linspace(-2.4, 2.4, num=25),
            "njets"  : range(12),
             "vmult" : [0,4,8,12,14,16,18,20,22,24,28,32,36],
            "met"   : [0,20,35,50,70,90,110,130,150,150,185,185,250],
            "nEle" : [0,1],
            "nPho" : [0,1],
            "hadTowOverEM": np.linspace(0, 0.01, num = 64),
            "dR": np.linspace(0, 2, num = 128)}

genTypes = ["genEle", "genMu", "genTau", "genPho", "genJet"]

recoTypes = ["recoPho", "recoEle"]

regions = ["barrel" , "endcap"]
regions2 = ["bb", "eb", "be"]

#eVars = ["pt", "eta" , "njets" , "met", "nEle", "nPho", "hadTowOverEM"]
eVars = ["dR", "hadTowOverEM"]

finalStates = ["ee", "eg", "ge", "gg"]

def init_hists(dType):
  hist_dict = OrderedDict()

  hist_dict["num_events"] = TH1F("num_events", "num_events", 1, 0, 1)

  for region in regions2:
    for finalState in finalStates:
      hName = finalState + "_" + region + "_matrix"
      hist = TH2I(hName, hName, 5, 0, 5, 5, 0, 5) #Lead, Trail
      for i in range(5):
        hist.GetXaxis().SetBinLabel(i+1, genTypes[i])
        hist.GetYaxis().SetBinLabel(i+1, genTypes[i])

      #hist.SetCanExtend(ROOT.TH1.kAllAxes)
      hist_dict[hName] = hist
      

  for recoType in recoTypes:
    for genType in genTypes:
      for eVar in eVars:
        for region in regions:
          h_bins = array("d", varBins[eVar])
          hName = recoType + "_" + region + "_noMatch_" + genType + "_" + eVar
          hist_dict[hName] = TH1F( hName, hName, len(h_bins)-1, h_bins)

  return hist_dict


def load_hists(dType, filename):
  f = TFile.Open(filename)

  if not f:
    print("File unable to be opened: " + filename)
    quit()

  hists = OrderedDict()

  for key in f.GetListOfKeys():
    h = key.ReadObj()
    h.SetDirectory(0) #Write to RAM
    hists[h.GetName()] = h

  f.Close()
  return hists


def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def pcdiff(a,b):
  return 100*np.abs(a-b)/(a+b)/2


def ratioPoissonError(A, B):
    A_err = np.sqrt(A)
    B_err = np.sqrt(B)
    
    ep1 = (A_err/B)**2
    ep2 = (A*B_err/B**2)**2
    
    return np.sqrt( ep1 + ep2 )


cuts_2016 = { "barrel" : {"loose": { "hadTowOverEM"       : None,
                                "sigmaIetaIeta"      : None,
                                "pfChargedIsoRhoCorr": None,
                                "pfNeutralIsoRhoCorr": None,
                                "pfGammaIsoRhoCorr"  : None},

                     "medium": {"hadTowOverEM"       : None,
                                "sigmaIetaIeta"      : None,
                                "pfChargedIsoRhoCorr": None,
                                "pfNeutralIsoRhoCorr": None,
                                "pfGammaIsoRhoCorr"  : None},

                      "tight": {"hadTowOverEM"       : 0.0269,
                                "sigmaIetaIeta"      : 0.00994,
                                "pfChargedIsoRhoCorr": 0.202,
                                "pfNeutralIsoRhoCorr": (0.264,0.0148, 0.000017),
                                "pfGammaIsoRhoCorr"  : (2.362, 0.0047)}},


         "endcap" : { "loose": { "hadTowOverEM"      : None,
                                "sigmaIetaIeta"      : None,
                                "pfChargedIsoRhoCorr": None,
                                "pfNeutralIsoRhoCorr": None,
                                "pfGammaIsoRhoCorr"  : None },

                     "medium": {"hadTowOverEM"       : None,
                                "sigmaIetaIeta"      : None,
                                "pfChargedIsoRhoCorr": None,
                                "pfNeutralIsoRhoCorr": None,
                                "pfGammaIsoRhoCorr"  : None },

                      "tight": {"hadTowOverEM"       : 0.0213,
                                "sigmaIetaIeta"      : 0.03,
                                "pfChargedIsoRhoCorr": 0.034,
                                "pfNeutralIsoRhoCorr": (0.586, 0.0163, 0.000014),
                                "pfGammaIsoRhoCorr"  : (2.617, 0.0034)} }}

def passID( era, ID, t, ipho):
 
  if not "16" in era:
    "Using wrong iso cuts for " + era
    quit()

  if t.Photons_isEB[ipho]:
    region = "barrel"
  else:
    region = "endcap"

  if t.Photons_hadTowOverEM[ipho] > cuts_2016[region][ID]["hadTowOverEM"]:
    return False

  if t.Photons_sigmaIetaIeta[ipho] > cuts_2016[region][ID]["sigmaIetaIeta"]:
    return False

  if t.Photons_pfChargedIsoRhoCorr[ipho] > cuts_2016[region][ID]["pfChargedIsoRhoCorr"]:
    return False

  pt = t.Photons[ipho].Pt()

  c = cuts_2016[region][ID]["pfNeutralIsoRhoCorr"]
  
  if t.Photons_pfNeutralIsoRhoCorr[ipho] > c[0] + c[1]*pt + c[2]*pt**2:
    return False

  c = cuts_2016[region][ID]["pfGammaIsoRhoCorr"]

  if t.Photons_pfGammaIsoRhoCorr[ipho] > c[0] + c[1]*pt:
    return False

  return True


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
