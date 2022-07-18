from __future__ import division
import numpy as np


nZ = 100000

#Define regions that we would see in data/mc
regions = ["BB", "EB", "BE", "EE"] ## First letter is the region of the electron, second is photon

selections = ["ee", "eg" , "gg"]

cDict = { region : { selection: 0 for selection in selections} for region in regions }

#Arbitrary probabilities
f_B = 0.015 ## Fakerate Barrel 
f_E = 0.025 ## Fakerate Endcap
eProb = 0.35 ## Probability of Object -> Endcap

class Electron:
  def __init__(self):
    self.barrel = np.random.rand() > eProb
    if self.barrel:
      self.xseed = np.random.rand() > f_B
    else:
      self.xseed = np.random.rand() > f_E

for event in range(nZ):
  lead = Electron()
  trail = Electron()

  if (lead.xseed and trail.xseed): #Passing probes ee

    if lead.barrel and trail.barrel:
      r = ["BB"]
    elif lead.barrel != trail.barrel:
      r = ["EB","BE"] #Not distinguishing between electrons
    else:
      r = ["EE"]

    for region in r:
      cDict[region]["ee"] += 1

  elif (lead.xseed != trail.xseed): #Failing probes eg
    electron = lead if lead.xseed else trail
    photon = trail if lead.xseed else lead

    region = ("B" if electron.barrel else "E") + ("B" if photon.barrel else "E")

    cDict[region]["eg"] += 1


#
print("The barrel fakerate is : " + str(f_B))
print("The endcap fakerate is : " + str(f_E))

print("The theoretical transfer factor is f/(1-f) Barrel: " + str(f_B/(1-f_B)))
print("The theoretical transfer factor is f/(1-f) Endcap: " + str(f_E/(1-f_E)))


for region in regions:
  print("region: " + region +" tf : " + str(cDict[region]["eg"]/cDict[region]["ee"]))

"""
print("Now generate EG events and try to predict GG from EG")

regions = ["B" , "E"]
selections = ["eg" , "gg"]
egDict = {region:{selection:0 for selection in selections} for region in regions}

nEG = 100000

for n in range(nEG):
  electron = Electron()
  
  region = "B" if electron.barrel else "E"
  selection = "eg" if electron.xseed else "gg"

  egDict[region][selection] += 1


print("The prediction for nGG is: nGG = nEG*tf")
for region in regions:
  print( region + " nEG = " + str(egDict[region]["eg"]))
  print( region + " nGG = " + str(egDict[region]["gg"]))
  tf = tf_barrel if region == "B" else tf_endcap
  print( region + " predicted GG = " + str(tf*egDict[region]["eg"]))
"""
