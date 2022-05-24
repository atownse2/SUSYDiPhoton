import numpy as np
import pickle
from ROOT import Math


def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)




def fakerate(eeSig, eeSigError, egSig, egSigError):
  fakerate = egSig/eeSig

  EPeg = 1/(eeSig)
  EPee = egSig/((eeSig**2))

  error = np.sqrt((EPeg*egSigError)**2 + (EPee*eeSigError)**2)

  return fakerate, error

def fakerate2(eeSig, eeSigError, egSig, egSigError):
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
