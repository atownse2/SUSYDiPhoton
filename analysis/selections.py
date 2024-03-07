import numpy as np

selections = {}

# Signal Region
selections["SR"] = {
    "PhotonWP" : "loose",
    "NPhotons" : 2,
    "PixelSeedVeto" : True,
    "MinBarrelPhotons" : 1,
    "NGoodJets" : 2,
    "HardMETPt" : 130,
    "AbsHardMetMinusMet" : 100,              
    "MuonVeto" : True,
    "ElectronVeto" : True,
}

# EGamma Control Region
selections["EG"] = selections["SR"].copy()
selections["EG"].pop("PixelSeedVeto")

# DY Control Region
selections["DY"] = {
    "PhotonWP" : "loose",
    "NPhotons" : 2,
    "ZWindow" : True,
}

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?sortcol=0;table=1;up=0#sorted_table
photonWPs = { 'loose' : {'barrel': { 'hadTowOverEM'       : 0.04596,
                                   'sigmaIetaIeta'      : 0.0106,
                                   'pfChargedIsoRhoCorr': 1.694,
                                   'pfNeutralIsoRhoCorr': (24.032, 0.01512, 2.259e-05),
                                   'pfGammaIsoRhoCorr'  : (2.876, 0.004017)}, 
                       'endcap': { 'hadTowOverEM'       : 0.0590,
                                   'sigmaIetaIeta'      : 0.0272,
                                   'pfChargedIsoRhoCorr': 2.089,
                                   'pfNeutralIsoRhoCorr': (19.722, 0.0117, 2.3e-05),
                                   'pfGammaIsoRhoCorr'  : (4.162, 0.0037)}}, 
            'medium': {'barrel': { 'hadTowOverEM'       : 0.02197,
                                   'sigmaIetaIeta'      : 0.01015,
                                   'pfChargedIsoRhoCorr': 1.141,
                                   'pfNeutralIsoRhoCorr': (1.189, 0.01512, 2.259e-05),
                                   'pfGammaIsoRhoCorr'  : (2.08, 0.004017)},
                        'endcap':{ 'hadTowOverEM'       : 0.0326,
                                   'sigmaIetaIeta'      : 0.0272,
                                   'pfChargedIsoRhoCorr': 1.051,
                                   'pfNeutralIsoRhoCorr': (2.718, 0.0117, 2.3e-5),
                                   'pfGammaIsoRhoCorr'  : (3.867, 0.0037)}}}

def pass_photon_ID(
    photonWP,
    pt,
    region,
    hadTowOverEM,
    sigmaIetaIeta,
    pfChargedIsoRhoCorr,
    pfNeutralIsoRhoCorr,
    pfGammaIsoRhoCorr):
    
    c = photonWPs[photonWP][region]
    if hadTowOverEM > c['hadTowOverEM']: return False
    if sigmaIetaIeta > c['sigmaIetaIeta']: return False
    if pfChargedIsoRhoCorr > c['pfChargedIsoRhoCorr']: return False

    p = c['pfNeutralIsoRhoCorr']
    if pfNeutralIsoRhoCorr > (p[0] + p[1]*pt + p[2]*pt*pt) : return False
    p = c['pfGammaIsoRhoCorr']
    if pfGammaIsoRhoCorr > (p[0] + p[1]*pt): return False

    return True

class EventSelection:

    def __init__(
        self,
        tree, 
        trigger_index,
        analysis_region,
        cutflow,
        ):

        cutflow['total'] += 1

        self.selections = selections[analysis_region]
        self.pass_selection = self.pass_selections(tree, trigger_index, cutflow)

    def pass_selections(self, tree, trigger_index, cutflow):
        s = self.selections
        # Preselections
        if tree.IsRandS != 0: return False
        if tree.TriggerPass[trigger_index]==0: return False
        cutflow['Trigger'] += 1

        # if tree.mva_Photons1Et<=75: return False # Don't know why this is here

        if "HardMETPt" in s:
            if tree.HardMETPt<=s["HardMETPt"]: return False
            cutflow['HardMETPt'] += 1

        if "NPhotons" in s:
            if len(tree.Photons) < s["NPhotons"]: return False
            cutflow["NPhotons"] += 1

        if 'NGoodJets' in s:
            if tree.mva_Ngoodjets<s['NGoodJets']: return False
            cutflow['NJets'] += 1
        
        candidates = self.get_photon_candidates(tree)

        if "NPhotons" in s:
            if len(candidates) < s["NPhotons"]: return
            cutflow['NLoosePhotons'] += 1

        # Analyis selections
        candidates = candidates[:2]
        if "MuonVeto" in s:
            if len(tree.Muons) > 0: return
            cutflow['MuonVeto'] += 1
        if "ElectronVeto" in s:
            hardElectrons = [e for e in tree.Electrons if e.Pt() > 40 and not any(deltaR(e, g['4vec']) < 0.1 for g in candidates)]
            if len(hardElectrons) > 0: return
            cutflow['ElectronVeto'] += 1
        
        if "MinBarrelPhotons" in s:
            NBarrelPho = sum([c['barrel'] for c in candidates])
            if NBarrelPho < s["MinBarrelPhotons"]: return
            cutflow['MinBarrelPhotons'] += 1

        if "ZWindow" in s:
            invMass = (candidates[0]['4vec'] + candidates[1]['4vec']).M()
            if invMass < 50 or invMass > 130: return
            cutflow['ZWindow'] += 1

        self.pass_selection = True
        self.candidates = candidates[:2]
        cutflow['All'] += 1

        return True

    def get_photon_candidates(
        self,
        tree,
        ):

        s = self.selections

        candidates = []
        for i, photon in enumerate(tree.Photons):
            candidate = {}

            if "PixelSeedVeto" in s and tree.Photons_hasPixelSeed[i]: continue
            if photon.Pt() < 80 or abs(photon.Eta()) > 2.4: continue

            id_vars = [
                photon.Pt(),
                'barrel' if tree.Photons_isEB[i] else 'endcap',
                tree.Photons_hadTowOverEM[i],
                tree.Photons_sigmaIetaIeta[i],
                tree.Photons_pfChargedIsoRhoCorr[i],
                tree.Photons_pfNeutralIsoRhoCorr[i],
                tree.Photons_pfGammaIsoRhoCorr[i],
            ]

            if pass_photon_ID('loose', *id_vars):
                candidate['photonWP'] = 'loose'
            elif pass_photon_ID('medium', *id_vars):
                candidate['photonWP'] = 'medium'
            else:
                continue

            candidate['index'] = i
            candidate['4vec'] = photon
            candidate['barrel'] = tree.Photons_isEB[i]
            candidate['xseed'] = bool(tree.Photons_hasPixelSeed[i])

            candidates.append(candidate)

        candidates = sorted( candidates, key = lambda i: i["4vec"].Pt(), reverse=True) 
        return candidates


def pcdiff(a,b):
    return 100*np.abs(a-b)/(a+b)

def deltaR(p1,p2):
    from ROOT import Math
    return Math.VectorUtil.DeltaR(p1,p2)

def genMatching(tree, reco_particles): 
    '''Match genParticles to EM objects'''
    for p in reco_particles:
        matches = []
        for iPar, genPar in enumerate(tree.GenParticles):

            dR = deltaR(p['4vec'], genPar)
            pcdiff_pt = pcdiff(p['4vec'].Pt(), genPar.Pt())
            if dR > 0.4 : continue

            # Classify gen particle
            pdgid = tree.GenParticles_PdgId[iPar]
            if isW(pdgid): continue
            elif isPhoton(pdgid) : genPar_name = 'genPho'
            elif isElectron(pdgid) : genPar_name = 'genEle'
            elif isMuon(pdgid) : genPar_name = 'genMu'
            elif isTau(pdgid) : genPar_name = 'genTau'
            else: genPar_name = 'genJet'

            matches.append({'index':iPar, 'name':genPar_name, 'pcdiff_pt':pcdiff_pt })

        if len(matches) > 1: # Sort matches by pt_pcdiff
            matches = sorted(matches, key = lambda i: i['pcdiff_pt'])
        else: # If no matches, call it a genJet
            matches.append({'index':None, 'name':'genJet', 'pcdiff_pt':None})

        # Add matches to reco particle container  
        p['genMatch'] = matches[0]

def isPhoton(pdgId):
    return abs(pdgId) == 22

def isNeutrino(pdgId):
    return abs(pdgId) in [12,14,16]

def isW(pdgId):
    return abs(pdgId) == 24

def isLepton(pdgId):
    return abs(pdgId) in [11,13,15]

def isElectron(pdgId):
    return abs(pdgId) == 11

def isMuon(pdgId):
    return abs(pdgId) == 13

def isTau(pdgId):
    return abs(pdgId) == 15

class WEvent:
    def __init__(self, tree):

        # Count W bosons
        self.nWs = 0
        for iPar, genPar in enumerate(tree.GenParticles):
            if isW(tree.GenParticles_PdgId[iPar]):
                self.nWs += 1
        
        self.hasW = self.nWs > 0
        if not self.hasW:
            return

        # Get daughters
        self.daughters = []
        for iPar, genPar in enumerate(tree.GenParticles):
            if isNeutrino(tree.GenParticles_PdgId[iPar]): continue
            iMom = tree.GenParticles_ParentIdx[iPar]
            if isW(tree.GenParticles_PdgId[iMom]):
                self.daughters.append({
                    'index': iPar,
                    '4vec': genPar,
                    'pdgId': tree.GenParticles_PdgId[iPar],
                })
  
        # Classify decay mode
        if any(isLepton(daughter['pdgId']) for daughter in self.daughters):
            self.decay_mode = 'leptonic'
        else:
            self.decay_mode = 'hadronic'

        # Check if daughters match with a pf candidate
        for daughter in self.daughters:
            daughter['hasMatch'] = False
            for jet in tree.JetsAUX:
                if deltaR(daughter['4vec'], jet) < 0.1:
                    daughter['hasMatch'] = True
                    break

        # Count daughters with matches
        self.nMatchedDaughters = sum([daughter['hasMatch'] for daughter in self.daughters])
        

if __name__ == '__main__':
  pass