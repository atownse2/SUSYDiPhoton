import sys
import os
import json
import re

from array import array
from collections import OrderedDict

from ROOT import Math
import numpy as np

top_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/'

# gROOT.SetBatch()        # don't pop up canvases
# gROOT.SetStyle('Plain') 



###Cuts
def make_cuts(tree, trigger_index, photonWP, skimmer):
  '''Reject events that don't pass cuts. Returns two highest pt photons/electrons which pass all cuts.'''
  if tree.TriggerPass[trigger_index]!=1 or tree.IsRandS != 0: return False
  if tree.mva_Ngoodjets<=1 or tree.mva_Photons1Et<=75: return False
  if skimmer == 'closure' and (tree.HardMETPt<=130 or abs(tree.HardMetMinusMet)>=100): return False
  if len(tree.Muons)   > 0: return False
  if len(tree.Photons) < 2: return False

  #Get photons that pass
  photonID = [passID(photonWP, tree.Photons_isEB[i], tree.Photons[i].Pt(),   
                tree.Photons_hadTowOverEM[i], 
                tree.Photons_sigmaIetaIeta[i], 
                tree.Photons_pfChargedIsoRhoCorr[i], 
                tree.Photons_pfNeutralIsoRhoCorr[i], 
                tree.Photons_pfGammaIsoRhoCorr[i]) for i in range(len(tree.Photons))]

  pho_pass = lambda i: photonID[i] and tree.Photons[i].Pt() > 80 and abs(tree.Photons[i].Eta()) < 2.4
  eg_pass = [{"pt"  : tree.Photons[i].Pt(),
              "eta" : tree.Photons[i].Eta(),
              "4vec"  : tree.Photons[i] , 
              "xseed" : bool(tree.Photons_hasPixelSeed[i]),
              "barrel": bool(tree.Photons_isEB[i]),
              "index" : i} 
              for i in range(len(tree.Photons)) if pho_pass(i)]

  if len(eg_pass) < 2:
    return False

  #Sort by pt
  eg = sorted( eg_pass, key = lambda i: i["pt"], reverse=True) 

  #At least one photon in the barrel
  if skimmer=='closure' and not eg[0]['barrel'] and not eg[1]['barrel']: return False

  #Veto if extra leptons
  #if skimmer=='closure' and len(tree.Electrons) > len([e for e in eg if e['xseed']]): return False 

  return eg[0],eg[1]

#WP cutBasedElectronID_Fall17_94X_V2
cutDict = { 'barrel' : {'loose': { 'hadTowOverEM'       : 0.04596,
                                   'sigmaIetaIeta'      : 0.0106,
                                   'pfChargedIsoRhoCorr': 1.694,
                                   'pfNeutralIsoRhoCorr': (24.032, 0.01512, 2.259e-05),
                                   'pfGammaIsoRhoCorr'  : (2.876, 0.004017)}, 
                        'medium':{ 'hadTowOverEM'       : 0.02197,
                                   'sigmaIetaIeta'      : 0.01015,
                                   'pfChargedIsoRhoCorr': 1.141,
                                   'pfNeutralIsoRhoCorr': (1.189, 0.01512, 2.259e-05),
                                   'pfGammaIsoRhoCorr'  : (2.08, 0.004017)}},
            'endcap' : {'loose': { 'hadTowOverEM'       : 0.0590,
                                   'sigmaIetaIeta'      : 0.0272,
                                   'pfChargedIsoRhoCorr': 2.089,
                                   'pfNeutralIsoRhoCorr': (19.722, 0.0117, 2.3e-05),
                                   'pfGammaIsoRhoCorr'  : (4.162, 0.0037)}, 
                        'medium':{ 'hadTowOverEM'       : 0.0326,
                                   'sigmaIetaIeta'      : 0.0272,
                                   'pfChargedIsoRhoCorr': 1.051,
                                   'pfNeutralIsoRhoCorr': (2.718, 0.0117, 2.3e-5),
                                   'pfGammaIsoRhoCorr'  : (3.867, 0.0037)}}}

def passID(WP, isEB, pt, hadTowOverEM, sigmaIetaIeta, pfChargedIsoRhoCorr, pfNeutralIsoRhoCorr, pfGammaIsoRhoCorr):
  region = 'barrel' if isEB else 'endcap'
  c = cutDict[region][WP]
  if hadTowOverEM > c['hadTowOverEM']:
    return False
  elif sigmaIetaIeta > c['sigmaIetaIeta']:
    return False
  elif pfChargedIsoRhoCorr > c['pfChargedIsoRhoCorr']:
    return False

  i = c['pfNeutralIsoRhoCorr']
  if pfNeutralIsoRhoCorr > (i[0] + i[1]*pt + i[2]*pt*pt) :
    return False

  i = c['pfGammaIsoRhoCorr']
  if pfGammaIsoRhoCorr > (i[0] + i[1]*pt):
    return False

  return True


def pcdiff(a,b):
  return 100*np.abs(a-b)/(a+b)/2

def deltaR(photon,jet):
  return Math.VectorUtil.DeltaR(photon, jet)

def genMatching(reco_g, genParticles, genParticles_PdgId): 
  '''Match genParticles to EM objects'''
  for eg in reco_g:
    hasMatch = False
    for iPar, genPar in enumerate(genParticles):
      if deltaR(eg['4vec'], genPar) < 1:
        pdgid = abs(genParticles_PdgId[iPar])

        if pdgid == 11:
          genPar_name = "genEle"
        elif pdgid == 13:
          genPar_name = "genMu"
        elif pdgid == 15:
          genPar_name = "genTau"
        elif pdgid == 22:
          genPar_name = "genPho"
        else:
          genPar_name = "genJet"

        #Match or not
        if hasMatch == False:
          eg['genPar_index'] = iPar
          eg['genPar_name'] = genPar_name
          eg['genPar_pt'] = genPar.Pt()
          hasMatch = True
        elif pcdiff(eg['pt'], genPar.Pt()) < pcdiff(eg['pt'], eg['genPar_pt']):
          eg['genPar_index'] = iPar
          eg['genPar_name'] = genPar_name
          eg['genPar_pt'] = genPar.Pt()

    if hasMatch == False:
      eg['genPar_name'] = 'genJet'
    
  return reco_g

def isNeutrino(pdgId):
  return abs(pdgId) in [12,14,16]

def isW(pdgId):
  return abs(pdgId) == 24

def isLepton(pdgId):
  return abs(pdgId) in [11,13,15]

def lostLepton(electrons, muons, genparticles, pdgids, genparticles_parentids):
  '''Check if there is a reconstructed electron or muon from a W decay'''
  recos = [y for x in [electrons, muons] for y in x]

  genLeptonsFromW = []
  for genPar, pdgid, genparticle_parentid in zip(genparticles, pdgids, genparticles_parentids):
    if isLepton(pdgid) and isW(genparticle_parentid):
      genLeptonsFromW.append(genPar)
  
  #If there are no leptons from W decay, return False
  for genLepton in genLeptonsFromW:
    hasRecoMatch  = False
    for reco in recos:
      if deltaR(genLepton, reco) < 0.6:
        hasRecoMatch = True
        break
    if not hasRecoMatch :
      return True
  return False


if __name__ == '__main__':
  pass