import os
import numpy as np
import XRootD.client
import subprocess
import re
import json

import sys
top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton"
sys.path.append(top_dir)

from analysis.logger import Logger

l = Logger()

eras = ["2016", "2017", "2018"]

mc_eras = {"2016" : "Summer16v3",
           "2017" : "Fall17",
           "2018" : "Autumn18"}

data_eras = {"2016" : "Run2016",
            "2017" : "Run2017",
            "2018" : "Run2018"}

era_lumi = {"2016" : 35.9, 
            "2017" : 41.5,
            "2018" : 59.6}

phases = { "Phase0": ["2016"],
           "Phase1": ["2017", "2018"]}

dType_tags_include_all = {'data' : ['Run201'],
                      'DYJetsToLL' : ['DYJetsToLL_M-50_HT'],
                      'WGJets' : ['WGJets_MonoPhoton_PtG'],
                      'TTGJets' : ['TTGJets'],}

dType_tags_include_any = {'data' : ["EGamma", "DoubleEG"],
                          'DYJetsToLL' : [""],
                      'WGJets' : [""],
                      'TTGJets' : [""],}

dType_tags_exclude_any = {'data' : [],
                      'DYJetsToLL' : [],
                      'WGJets' : [],
                      'TTGJets' : [],}

dTypes = list(dType_tags_include_any.keys())

sample_subsets = {"DYJetsToLL" : ["DYJetsToLL_M-50_HT-100to200", 
                                  "DYJetsToLL_M-50_HT-200to400",
                                  "DYJetsToLL_M-50_HT-400to600",
                                  "DYJetsToLL_M-50_HT-600to800",
                                  "DYJetsToLL_M-50_HT-800to1200",
                                  "DYJetsToLL_M-50_HT-1200to2500",
                                  "DYJetsToLL_M-50_HT-2500toInf"],
                  "WGJets": ["WGJets_MonoPhoton_PtG-40to130", 
                             "WGJets_MonoPhoton_PtG-130"],
                  "TTGJets": ["TTGJets"]}

trigger_names = {'2016' : 'HLT_DoublePhoton60_v',
                  '2017' : 'HLT_DoublePhoton70_v',
                  '2018' : 'HLT_DoublePhoton70_v'}

top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton"
skim_dir = "/store/group/lpcsusyphotons/"
filelist_dir = f"{top_dir}/analysis/metadata/filelists"

def tag_filename(dType, filename):
  if all([tag in filename for tag in dType_tags_include_all[dType]]) and \
     any([tag in filename for tag in dType_tags_include_any[dType]]) and \
      not any([tag in filename for tag in dType_tags_exclude_any[dType]]):
    return True
  else:
    return False

def get_dType(filename):
  for dType in dTypes:
    if tag_filename(dType, filename):
      return dType
  
  raise ValueError('dType not found for {}'.format(filename))

def get_subset(filename):
  dType = get_dType(filename)
  for subset in sample_subsets[dType]:
    if subset in filename:
      return subset
  raise ValueError('sample subset not found for {}'.format(filename))

def get_era(filename):
  dType = get_dType(filename)
  if dType == 'data':
    for era, era_tag in data_eras.items():
      if era_tag in filename:
        return era
  else:
    for era, era_tag in mc_eras.items():
      if era_tag in filename:
        return era
  raise ValueError('era not found for {}'.format(filename))

### XrootD functions
def xrd_dirlist(path, redirector="ndcms.crc.nd.edu"):
  # Get list of directories using XrootD
  dirlist = []
  fs = XRootD.client.FileSystem(redirector)
  status, listing = fs.dirlist(path)
  print(status)
  for f in listing:
    dirlist.append(f.name)
  return dirlist

# Get and cache the file lists
def get_filelist(skim, remake=False):
  filelist_name = 'filelists/{}_filelist.txt'.format(skim)
  if not os.path.exists(filelist_name) and not remake:
    # Use XrootD to get list of files
    filelist = xrd_dirlist(skim_dir + skim)
    # Write to file
    with open(filelist_name, 'w') as f:
      for filename in filelist:
        f.write(filename+'\n')
  else:
    # Read from file
    with open(filelist_name, 'r') as f:
      filelist = [f.replace('\n','') for f in f.readlines()]
  return filelist
###


def get_ntuple_filenames(dType, skim, 
                         runbatch=False, batch=None, batches=None, 
                         redirector="ndcms.crc.nd.edu",
                         test=False):
  '''Get input filenames'''
  #redirector = "deepthought.crc.nd.edu"
  #redirector = "cmsxrootd.fnal.gov"
  ntuple_location = f"root://{redirector}//store/group/lpcsusyphotons/{skim}"
  txtfile = f"{filelist_dir}/{skim}_filelist.txt"

  filenames = []
  with open(txtfile, "r") as f:
    for line in f:
      if tag_filename(dType, line):
        filenames.append(f"{ntuple_location}/{line[:-1]}")

  if test:
    # Write filenames to tmp file in cache
    tmp_file = f"{top_dir}/test/tmp/{dType}_{skim}_filelist.txt"
    with open(tmp_file, "w") as f:
      f.write("\n".join(filenames))

  if runbatch:
    # Assert that batch and batches are defined
    assert batch is not None
    assert batches is not None

    file_batches = np.array_split(filenames, batches)
    filenames = file_batches[batch]

  if not any(filenames):
    print("No files in batch")
    quit()
  
  return filenames


def characterize_skims():
  # Get the number of files for each dType and era in each skim and store it in a pandas dataframe with columns: dType_era and rows: skim
  import pandas as pd
  
  tags = ["Skims", "TreeMaker","skims"]
  skims = [ d for d in xrd_dirlist(skim_dir) if any(tag in d for tag in tags) ]
  df = pd.DataFrame(index=skims, columns=['{}_{}'.format(dType, era) for dType in dTypes for era in eras])
  for skim in skims:
    filelist = get_filelist(skim)
    for dType in dTypes:
      for era in eras:
        df.loc[skim, '{}_{}'.format(dType, era)] = len([f for f in filelist if dType in f and era in f])

  # Save dataframe to YAML
  import yaml

  # Convert the dataframe to a dictionary
  data = df.transpose().to_dict()

  # Save the dictionary to a YAML file
  with open('skim_info.yaml', 'w') as f:
      yaml.dump(data, f, default_flow_style=False)


def get_trigger_index(filename, trigger_name, tree_name='TreeMaker2/PreSelection', batch=0):

  dType = get_dType(filename)
  era = get_era(filename)

  skims = [ f.replace("_filelist.txt","") for f in os.listdir(filelist_dir) if '_filelist.txt' in f]
  skim = [s for s in skims if s in filename][0] # Should be only one match anyway

  trigger_dict = {}
  cached = False
  trigger_file = f'{top_dir}/analysis/metadata/cache/trigger_index.json'
  if os.path.exists(trigger_file):
    with open(trigger_file, 'r') as f:
      trigger_dict = json.load(f)
      if (dType in trigger_dict.keys()) and \
          (skim in trigger_dict[dType].keys()) and \
          (era in trigger_dict[dType][skim].keys()):
        cached = True

  if cached:
    trigger_index = trigger_dict[dType][skim][era]
  else:
    # Print trigger names
    output = get_tree_print(filename, tree_name, batch)
    
    # Parse the output
    trigger_names = re.findall(r'HLT_[\w\d_]*', output.decode('utf-8'))
    trigger_index = trigger_names.index(trigger_name)

    # Save to cache
    if dType not in trigger_dict.keys():
      trigger_dict[dType] = {}
    if skim not in trigger_dict[dType].keys():
      trigger_dict[dType][skim] = {}
    trigger_dict[dType][skim][era] = trigger_index
    with open(trigger_file, 'w') as f:
      json.dump(trigger_dict, f, indent=2)
  
  l.log('Trigger index for {} in {} and era {} is {}'.format(trigger_name, skim, era, trigger_index),1)
  return trigger_index

def cache_trigger_indicies(skim, dTypes):
  for dType in dTypes:
    filenames = get_ntuple_filenames(dType, skim)
    isMC = dType != 'data'
    # Sort by era into dict
    era_tags = mc_eras if isMC else data_eras
    era_dict = {}
    for era in eras:
      for filename in filenames:
        if era_tags[era] in filename:
          era_dict[era] = filename
          break

    # Get trigger index for each era
    for era, filename in era_dict.items():
        trigger_name = trigger_names[era]
        get_trigger_index(filename, trigger_name)


def get_tree_print(filename, tree_name, batch):
    # Write a temporary ROOT macro
    macro_dir = '/tmp/atownse2/'
    macro_name = 'tmp_{}'.format(batch)
    macro_path = macro_dir + macro_name + '.C'

    with open(macro_path, 'w') as f:
        f.write('''
                void %s() {
                    TFile *f = TFile::Open("%s");
                    TTree *tree = (TTree*)f->Get("%s");
                    TBranch *branch = tree->GetBranch("TriggerPass");
                    branch->Print();
                    f->Close();
                }
                '''% (macro_name, filename, tree_name))
  
    # Run the macro with ROOT and capture the output
    output = subprocess.check_output(['root', '-l', '-b', '-q', macro_path])
    # Clean up
    subprocess.call(['rm', macro_path])
    return output

if __name__ == '__main__':
  l.set_max_verbosity(2)
  #characterize_skims()
  file = 'root://ndcms.crc.nd.edu//store/group/lpcsusyphotons/Skims_loosePhotonV10FullId/posterior-Run2016B-17Jul2018_ver2-v1.DoubleEG_0_RA2AnalysisTree_part1.root'
  #file = 'root://ndcms.crc.nd.edu//store/group/lpcsusyphotons/Skims_mediumPhotonV10/posterior-Run2016B-17Jul2018_ver2-v1.DoubleEG_0_RA2AnalysisTree_part1.root'

  #get_trigger_index(file, 'HLT_DoublePhoton60_v')
  cache_trigger_indicies('skimsforAustin', ["data", "DYJetsToLL"])
  cache_trigger_indicies('Skims_mediumPhotonV10', ["data", "WGJets", "TTGJets"])