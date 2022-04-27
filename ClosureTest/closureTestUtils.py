from ROOT import TH1F
from array import array
import sys
import numpy as np


#barrel tagHpt fakerate: 0.0163725060981 error: 0.00010434581740764053
#barrel tagLpt fakerate: 0.0154553105352 error: 0.0001780002362461558
#endcap tagHpt fakerate: 0.0229211070391 error: 0.0036137292754381896
#endcap tagLpt fakerate: 0.0205385071674 error: 0.003527769324111087



fakerateDict = {"Summer16v3" : { "barrel" : {"tagHpt" : (20.0e-03/2, 10e-05),
                                             "tagLpt" : (20.0e-03/2, 17e-05)},
                                 "endcap" : {"tagHpt" : (21.0e-03/2, 360e-05),
                                             "tagLpt" : (21.0e-03/2, 352e-05)} }}


def getHistName(sel, hType, ptRange, region = False, ptTag = False, bdt = False):
  hName = sel + "_" + hType + "_"

  if bdt:
    hName += bdt + "_"

  if region:
    hName += region + "_"

  if ptTag:
    hName += ptTag + "_"

  return hName + ptRange

### Define hist dicts and bins
ptRangeDict = {"WGJets" :["PtG-40to130", "PtG-130"] ,
               "DYJetsToLL_M-50": ["100to200","200to400","400to600","600to800","800to1200","1200to2500","2500toInf"],
               "TTJets": ["SingleLeptFromTbar", "SingleLeptFromT", "DiLept"]}

varBins = { "pt"    : array("d", np.round(np.linspace(80,300, num = 32), decimals = 3)),
            "eta"   : array("d", np.round(np.linspace(0,2.4, num = 32), decimals = 2)),
            "met"   : array("d", [0,10,20,35,50,70,90,110,130,150,150,185,185,250]),
            "DiEMpt": array("d", np.round(np.linspace(0, 500, num = 64), decimals = 1))}

def drawOverflow(h):
  h.Sumw2()

  nx = h.GetNbinsX() + 1
  xbins = array("d", [0]*(nx+1))
  for i in range(nx):
    xbins[i] = h.GetBinLowEdge(i+1)
  xbins[nx] = xbins[nx-1] + h.GetBinWidth(nx)
  tmpName = h.GetName() + "Overflow"
  htmp = TH1F(tmpName, h.GetName(), nx, xbins)

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

