from utils import pcdiff, deltaR, getHistName

def make_cuts(t, hists, HLT):
  if t.TriggerPass[HLT] != 1:
    return 0

  if len(t.Photons) < 2: 
    return 0

  EMpass = [ {"pt"  : t.Photons[i].Pt(),
              "eta" : t.Photons[i].Eta(),
              "4vec"  : t.Photons[i] , 
              "xseed" : bool(t.Photons_hasPixelSeed[i]),
              "barrel": t.Photons_isEB[i],
              "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 80 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

  if len(EMpass) < 2:
    return 0

  em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first

  if not t.Photons_isEB[em[0]["index"]]  and not t.Photons_isEB[em[1]["index"]]:
    return 0
      
  clean_jets = []

  for jet in [j for j in t.Jets if j.Pt() > 30]:
    jetPass = True
    for e in em:
      if deltaR(jet, e["4vec"]) < 0.4:
        jetPass = False
        break
    if jetPass:
      clean_jets.append(jet)

  njets = len(clean_jets)

  vmult = t.NVtx

  #w = t.CrossSection*t.puWeight*1000*35.9
  w = 1

  for reco in em:
    if reco["barrel"]:
      region = "barrel"
    else:
      region = "endcap"

    for iPar, genPar in enumerate(t.GenParticles):
      if genPar.Pt() < 10:
        continue

      if deltaR(reco["4vec"], genPar) > 0.2 and pcdiff(reco["pt"], genPar.Pt()) > 20:
        continue

      if abs(t.GenParticles_PdgId[iPar]) == 11:
        genPar_name = "genEle"  
      elif abs(t.GenParticles_PdgId[iPar]) == 15:
        genPar_name = "genTau"
      elif abs(t.GenParticles_PdgId[iPar]) == 22:
        genPar_name = "genPho"
      else:
        genPar_name = "genJet"

      if reco["xseed"]: #Electron
        #Fill Hists
        hists[getHistName(genPar_name, "recoEle", region, "pt" )].Fill(reco["4vec"].Pt() , w)
        hists[getHistName(genPar_name, "recoEle", region, "eta")].Fill(reco["4vec"].Eta(), w)
        hists[getHistName(genPar_name, "recoEle", region, "njets")].Fill(njets, w)
        hists[getHistName(genPar_name, "recoEle", region, "met")].Fill(t.MET, w)          
        hists[getHistName(genPar_name, "recoEle", region, "nEle")].Fill(0, w)

      else: #Reco Photon
        #Fill Hists
        hists[getHistName(genPar_name, "recoPho", region, "pt" )].Fill(reco["4vec"].Pt() , w)
        hists[getHistName(genPar_name, "recoPho", region, "eta")].Fill(reco["4vec"].Eta(), w)
        hists[getHistName(genPar_name, "recoPho", region, "njets")].Fill(njets, w)
        hists[getHistName(genPar_name, "recoPho", region, "met")].Fill(t.MET, w)
        hists[getHistName(genPar_name, "recoPho", region, "nPho")].Fill(0, w)



def make_samscuts(t, hists, HLT):
  if t.TriggerPass[HLT] != 1:
    return    

  EMpass = [ {"pt"  : t.Photons[i].Pt(),
              "eta" : t.Photons[i].Eta(),
              "4vec"  : t.Photons[i] , 
              "xseed" : bool(t.Photons_hasPixelSeed[i]),
              "barrel": t.Photons_isEB[i],
              "index" : i} for i in range(len(t.Photons)) if t.Photons[i].Pt() > 30 and abs(t.Photons[i].Eta()) < 2.4 and t.Photons_fullID[i]]

  em = sorted( EMpass, key = lambda i: i["pt"], reverse=True) #Highest Pt Objects are first

  clean_jets = []

  for jet in [j for j in t.Jets if j.Pt() > 30]:
    jetPass = True
    for e in em:
      if deltaR(jet, e["4vec"]) < 0.4:
        jetPass = False
        break
    if jetPass:
      clean_jets.append(jet)

  njets = len(clean_jets)

  vmult = t.NVtx

  #w = t.CrossSection*t.puWeight*1000*35.9
  w = 1


  for reco in em:
    if reco["barrel"]:
      region = "barrel"
    else:
      region = "endcap"

    for iPar, genPar in enumerate(t.GenParticles):
      if genPar.Pt() < 10:
        continue

      if deltaR(reco["4vec"], genPar) > 0.2 and pcdiff(reco["pt"], genPar.Pt()) > 20:
        continue

      if abs(t.GenParticles_PdgId[iPar]) == 11:
        genPar_name = "genEle"  
      elif abs(t.GenParticles_PdgId[iPar]) == 15:
        genPar_name = "genTau"
      elif abs(t.GenParticles_PdgId[iPar]) == 22:
        genPar_name = "genPho"
      else:
        genPar_name = "genJet"

      if reco["xseed"]: #Electron
        #Fill Hists
        hists[getHistName(genPar_name, "recoEle", region, "pt" )].Fill(reco["4vec"].Pt() , w)
        hists[getHistName(genPar_name, "recoEle", region, "eta")].Fill(reco["4vec"].Eta(), w)
        hists[getHistName(genPar_name, "recoEle", region, "njets")].Fill(njets, w)
        hists[getHistName(genPar_name, "recoEle", region, "met")].Fill(t.MET, w)          
        hists[getHistName(genPar_name, "recoEle", region, "nEle")].Fill(0, w)

      else: #Reco Photon
        #Fill Hists
        hists[getHistName(genPar_name, "recoPho", region, "pt" )].Fill(reco["4vec"].Pt() , w)
        hists[getHistName(genPar_name, "recoPho", region, "eta")].Fill(reco["4vec"].Eta(), w)
        hists[getHistName(genPar_name, "recoPho", region, "njets")].Fill(njets, w)
        hists[getHistName(genPar_name, "recoPho", region, "met")].Fill(t.MET, w)
        hists[getHistName(genPar_name, "recoPho", region, "nPho")].Fill(0, w)






version_mapping = { "default" : make_cuts,
                    "samscuts": make_samscuts}
