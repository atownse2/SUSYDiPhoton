import os

from array import array
from collections import OrderedDict

from ROOT import TFile, TH1F, TH2F, gROOT
import numpy as np

from analysis.metadata import sample_info as s
from analysis import logger

import hist

l = logger.Logger()

# Define hist axes
mc_dataset_axis = { 'WGJets': hist.axis.StrCategory([], name='mc_dataset'),}


gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') 

hists_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/hists'

final_states = ['eBeB', 'eBeE', 'eEeB', 'eEeE', # ee
                'gBgB', 'gBgE', 'gEgB', 'gEgE', # gg
                'eBgB', 'gBeB', #eg BB
                'eEgE', 'gEeE', #eg EE
                'eBgE', 'gEeB', 'eEgB', 'gBeE'] # eg EB

gg_final_states = [fs for fs in final_states if fs.count('g') == 2]
eg_final_states = [fs for fs in final_states if fs.count('g') == 1 and fs.count('e') == 1]


det_regions = ['barrel', 'endcap']



final_state_tf_region_map = { 'eBgB' : 'barrel',
                              'gBeB' : 'barrel',
                              'eBgE' : 'barrel',
                              'gBeE' : 'endcap',
                              'eEgB' : 'endcap',
                              'gEeB' : 'barrel'}

prediction_map = { 'gBgB' : ['eBgB', 'gBeB'],
                   'gBgE' : ['eBgE', 'gBeE'], 
                   'gEgB' : ['eEgB', 'gEeB']}

gg_regions = ['BB', 'BE', 'EB']

###


genEvents = ['genGJet', 'genEG']
genTypes = ["genEle", "genMu", "genTau", "genPho", "genJet"]
gg_cutflow = ['all', 'fakehardMET', 'realhardMET', 'realhardMET+lostLepton']

analysis_bins = {'BDT' : array("d",[-1, -0.13, 0.03, 1]),
                 'hardMET' : array("d",[130, 150, 185, 250]),
}

var_bins = {'bdt' : array("d", np.linspace(-1,1, num=32)),
        'pt' : array("d", [80,100,150,200,300]),
        'eta': array("d", [0, 0.8, 1.4442, 1.566, 2.0, 2.5]),
        'met' : array("d", [130, 150, 185, 250]),
        'inv mass' : array("d", np.linspace(60,120, num=32)),
}

var_bin_names = {'bdt' : ['lowBDT', 'medBDT', 'highBDT'],
              'pt' : ['pt0to{}'.format(int(var_bins['pt'][0]))] + \
                     ['pt{}to{}'.format(int(var_bins['pt'][i]),int(var_bins['pt'][i+1])) for i in range(len(var_bins['pt'])-1)] + \
                     ['pt{}toInf'.format(int(var_bins['pt'][-1]))],
              'met': ['hardMET0to{}'.format(int(var_bins['met'][0]))] + \
                     ['hardMET{}to{}'.format(int(var_bins['met'][i]),int(var_bins['met'][i+1])) for i in range(len(var_bins['met'])-1)] + \
                     ['hardMET{}toInf'.format(int(var_bins['met'][-1]))],
}

###Get bins
def get_BDT_bin(bdtVal):
  if bdtVal < -0.13:
    return "lowBDT"
  elif bdtVal < 0.03:
    return "medBDT"
  else:
    return "highBDT"
  

  ###Histogram
def write_hists(hists, dType, hist_type, version, runbatch=False, batch=None):
  outpath = get_hists_filename(dType, hist_type, version, runbatch=runbatch, batch=batch)
  print( "Saving hists to: " +  outpath)
  outfile = TFile.Open(outpath, 'RECREATE')
  if dType != 'data':
    for subset, histDict in hists.items():
      for histName, hist in histDict.items():
        hist.SetName(f"{subset}_{histName}")
        hist.Write()
  else:
    for histName, hist in hists.items():
      hist.Write()

  outfile.Close()


def get_hists_in_subsets(dType, hists):
  hists_in_subsets = {}
  for subset in s.sample_subsets[dType]:
    hists_in_subsets[subset] = {}
    for hName, h in hists.items():
      hists_in_subsets[subset][hName] = h.Clone()
  return hists_in_subsets


def load_hists(dType, skim_type, skim_version, verbosity=0):
  '''
  Loads histograms from root file into a dictionary. 
  If it is MC then scale by the total number of events.
  '''
  l.set_max_verbosity(verbosity)
  l.log(f"Loading {dType} hists...",1)
  filename = get_hists_filename(dType, skim_type, skim_version)
  if not os.path.exists(filename):
    l.log("Looking for batch files...",1)
    add_batch_hists(dType, skim_type, skim_version)

  isMC = dType != 'data'
  
  f = TFile.Open(filename)

  unscaled_hists = {}
  for key in f.GetListOfKeys():
    h = key.ReadObj()
    h.SetDirectory(0) #Write to RAM
    unscaled_hists[h.GetName()] = h

  if isMC:
    l.log(f"Scaling {dType} hists...",1)
    eras = s.eras
    if dType == "DYJetsToLL" and skim_type == "egtf":
      eras = eras[:2]
    mc_subsets = s.sample_subsets[dType]

    hists = {}
    for era in eras:
      l.log(f"Scaling {era}...",1)
      lumi = s.era_lumi[era]
      subset_0 = mc_subsets[0]
      n_events = unscaled_hists[f"{subset_0}_hEvents_{era}"].GetBinContent(1)
      l.log(f"n_events: {n_events}",2)
      scale = (lumi*1000)/n_events

      for hName, h in unscaled_hists.items():
        if subset_0 in hName and era in hName:
          hName = hName.replace(f"{subset_0}_","")
          h.SetName(hName)
          hists[hName] = h
          hists[hName].Scale(scale)
      
      if len(mc_subsets) > 1:
        for subset in mc_subsets[1:]:
          n_events = unscaled_hists[f"{subset}_hEvents_{era}"].GetBinContent(1)
          scale = (lumi*1000)/n_events
        for hName, h in unscaled_hists.items():
          if subset in hName and era in hName:
            hName = hName.replace(f"{subset}_","")
            hists[hName].Add(h, scale)

  else:
    hists = unscaled_hists

  f.Close()
  return hists

def add_batch_hists(dType, hist_type, version):
  filepath = get_hists_filename(dType, hist_type, version)
  if os.path.exists(filepath):
    print("Output file already exists: " + filepath)
    return

  filetag = get_hists_filename(dType, hist_type, version, runbatch=True, batch="").replace('.root','')
  batch_dir = os.path.dirname(filetag)

  if not any([filetag.split("/")[-1] in f for f in os.listdir(batch_dir)]):
    print("No batch files found for: " + filetag)
    raise Exception("No batch files found for: " + filetag)

  os.system(f'hadd -f {filepath} {filetag}*.root')
  #os.system(f'rm {filetag}*')

def get_hists_filename(dType, hist_type, version,
                       runbatch=False, batch=None):

  outfilename = f'{dType}_Run2_{hist_type}_{version}'
  outdir = hists_dir
  if runbatch:
    outfilename += "_{}".format(batch)
    outdir += '/batch'
  filepath = f"{outdir}/{outfilename}.root"

  return filepath


def background_matrix(hname):
  hist = TH2F(hname, hname, 5, 0, 5, 5, 0, 5) #Lead, Trail
  for i in range(5):
    hist.GetXaxis().SetBinLabel(i+1, genTypes[i])
    hist.GetYaxis().SetBinLabel(i+1, genTypes[i])
  return hist