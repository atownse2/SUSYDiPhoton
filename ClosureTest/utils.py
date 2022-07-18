from ROOT import TH1F, Math
from array import array
import sys
import numpy as np



fakerateDict = {'Summer16v3': {'barrel': {'tagHpt': (0.017320856027164676, 1.9904081469657285e-05), 'tagLpt': (0.017320856027164676, 1.9904081469657285e-05)}, 'endcap': {'tagHpt': (0.02502899833783366, 0.0014649742137513352), 'tagLpt': (0.02502899833783366, 0.0014649742137513352)}}}

def get_file_name(dType, era, ntuple_version, script_version, git_version, nbatch=None):
  f = dType + "_" + era + "_" + ntuple_version + "_" + script_version + "_" + git_version

  if nbatch is not None:
    f += "_" + str(nbatch)

  return f + ".root" 

def drawOverflow(h):
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

def getBDTbin(bdtVal):
  if bdtVal < -0.13:
    return "lowBDT"
  elif bdtVal < 0.03:
    return "medBDT"
  else:
    return "highBDT"


def deltaR(photon,jet):
  return Math.VectorUtil.DeltaR(photon, jet)
