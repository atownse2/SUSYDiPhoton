import os
import numpy as np

import re
import json

import sys
top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(top_dir)

from analysis.utils.logger import Logger

l = Logger()

years = ["2016", "2017", "2018"]

mc_years = {
  "2016" : "Summer16v3",
  "2017" : "Fall17",
  "2018" : "Autumn18"
  }

data_years = {
  "2016" : "Run2016",
  "2017" : "Run2017",
  "2018" : "Run2018"
  }

lumis = {
  "2016" : 35.9, 
  "2017" : 41.5,
  "2018" : 59.6
  }

phases = {
  "Phase0": ["2016"],
  "Phase1": ["2017", "2018"]
  }

dType_tags_include_all = {
  'data' : ['Run201'],
  'signal': ['SMS'],
  'DYJetsToLL' : ['DYJetsToLL_M-50_HT'],
  'WGJets' : ['WGJets_MonoPhoton_PtG'],
  'TTGJets' : ['TTGJets'],
  }

dType_tags_include_any = {
  'data' : ["EGamma", "DoubleEG"],
  'signal': ['T5Wg', 'T6Wg'],
  'DYJetsToLL' : [""],
  'WGJets' : [""],
  'TTGJets' : [""],
  }

dType_tags_exclude_any = {'data' : [],
                      'signal': [],
                      'DYJetsToLL' : [],
                      'WGJets' : [],
                      'TTGJets' : [],}

all_dTypes = ['data', 'signal', 'DYJetsToLL', 'WGJets', 'TTGJets']

sample_subsets = {
  "DYJetsToLL" : ["DYJetsToLL_M-50_HT-100to200", 
                  "DYJetsToLL_M-50_HT-200to400",
                  "DYJetsToLL_M-50_HT-400to600",
                  "DYJetsToLL_M-50_HT-600to800",
                  "DYJetsToLL_M-50_HT-800to1200",
                  "DYJetsToLL_M-50_HT-1200to2500",
                  "DYJetsToLL_M-50_HT-2500toInf"],
  "WGJets": ["WGJets_MonoPhoton_PtG-40to130", 
             "WGJets_MonoPhoton_PtG-130"],
  "TTGJets": "TTGJets",
  "data": "data",
  "signal": ["SMS-T5Wg", "SMS-T6Wg"],}

trigger_names = {
  "2016" : 'HLT_DoublePhoton60_v',
  "2017" : 'HLT_DoublePhoton70_v',
  "2018" : 'HLT_DoublePhoton70_v'
  }

top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton"
cache_dir = f"{top_dir}/cache"

filelist_dir = f"{cache_dir}/filelists"
skim_dir = "/store/group/lpcsusyphotons/"

def get_years_from_era(era):
  if era == "Run2":
    return ["2016", "2017", "2018"]
  elif era in years:
    return [era]
  elif era in phases:
    return phases[era]
  else:
    raise ValueError(f"Invalid era {era}")

def get_datasets(dType, year):
  if isinstance(sample_subsets[dType], list):
    return [f"{dType}_{subset}_{year}" for subset in sample_subsets[dType]]
  else:
    return [f"{dType}_{year}"]

class DataFile:
  def __init__(self, filename):
    self.filename = filename
    self.dType = get_dType(filename)
    self.skim_version = get_skim_version(filename)
    self.isMC = False if self.dType=='data' else True
    self.year = get_year(filename)
    self.dataset = get_subset(filename)
    self.trigger_index = get_trigger_index(self.dType, self.year, self.skim_version)

def tag_filename(dType, filename, years=years):

  if dType == 'data':
    year_tags = [data_years[year] for year in years]
  else:
    year_tags = [mc_years[year] for year in years]

  if not any([tag in filename for tag in year_tags]):
    return False

  if all([tag in filename for tag in dType_tags_include_all[dType]]) and \
     any([tag in filename for tag in dType_tags_include_any[dType]]) and \
      not any([tag in filename for tag in dType_tags_exclude_any[dType]]):
    return True
  else:
    return False

def get_dType(filename):
  for dType in all_dTypes:
    if tag_filename(dType, filename):
      return dType
  raise ValueError('dType not found for {}'.format(filename))

def get_subset(filename):
  dType = get_dType(filename)
  if dType == 'signal':
    match = re.search(r'(T[56]Wg)_m(\d+\.0)d(\d+\.0)', filename)
    if match:
      return "{}_{}_{}".format(match.group(1), int(float(match.group(2))), int(float(match.group(3))))
    else:
      raise ValueError('subset not found for {}'.format(filename))

  if isinstance(sample_subsets[dType], str): return sample_subsets[dType]
  for subset in sample_subsets[dType]:
    if subset in filename:
      return subset
  raise ValueError('subset not found for {}'.format(filename))

def get_year(filename):
  dType = get_dType(filename)
  if dType == 'data':
    for year, year_tag in data_years.items():
      if year_tag in filename:
        return year
  else:
    for year, year_tag in mc_years.items():
      if year_tag in filename:
        return year
  raise ValueError('year not found for {}'.format(filename))

def get_skim_version(filename):
  assert(skim_dir in filename)
  return filename.split('/')[-2]

def cache_filedict(ntuple_version, file_cache,
                   redirector="ndcms.crc.nd.edu"):

  ntuple_location = f"root://{redirector}//store/group/lpcsusyphotons/{ntuple_version}"
  txtfile = f"{filelist_dir}/{ntuple_version}_filelist.txt"

  if not os.path.exists(file_cache):
    filenames = [f"{ntuple_location}/{line[:-1]}" for line in open(txtfile, "r")]
    filedict = {}
    for filename in filenames:
      subset = get_subset(filename)
      if subset is None:
        continue
      year = get_year(filename)

      if subset not in filedict.keys():
        filedict[subset] = {}
      
      if year not in filedict[subset].keys():
        filedict[subset][year] = []
      
      filedict[subset][year].append(filename)

    with open(file_cache, "w") as f:
      json.dump(filedict, f, indent=2)

def get_filelist_from_filedict(dataset, year, ntuple_version,
                               nbatches=1, batch=0,
                               remove_bad_files=True):
  file_cache = f"{cache_dir}/{ntuple_version}_filedict.json"

  if not os.path.exists(file_cache):
    cache_filedict(ntuple_version, file_cache)
  
  with open(file_cache, "r") as f:
    filedict = json.load(f)
  
  if dataset not in filedict.keys():
    raise ValueError(f"Dataset {dataset} not found in filedict")
  if year not in filedict[dataset].keys():
    raise ValueError(f"Year {year} not found in filedict for dataset {dataset}")
  
  filenames = filedict[dataset][year]

  # if remove_bad_files:
  #   bad_filelist = f"{cache_dir}/bad_ntuples.json"
  #   if os.path.exists(bad_filelist):
  #     with open(bad_filelist, "r") as f:
  #       bad_file_dict = json.load(f)
  #     for bad_file in bad_file_dict['bad files']:
  #       for filename in filenames:
  #         if bad_file in filename:
  #           filenames.remove(filename)

  file_batches = np.array_split(filenames, nbatches)
  return file_batches[batch]


### XrootD functions
def xrd_dirlist(path, redirector="ndcms.crc.nd.edu"): # "cmsxrootd.fnal.gov"
  # Get list of directories using XrootD
  import XRootD.client
  dirlist = []
  fs = XRootD.client.FileSystem(redirector)
  status, listing = fs.dirlist(path)
  print(status)
  for f in listing:
    dirlist.append(f.name)
  return dirlist

# Get and cache the file lists
def get_filelist(skim_version, remake=False, max_verbosity=0):
  l.set_max_verbosity(max_verbosity)

  l.log(f"Getting filelist for {skim_version}", 1)
  filelist_name = f'{filelist_dir}/{skim_version}_filelist.txt'
  if not os.path.exists(filelist_name) and not remake:
    # Use XrootD to get list of files
    filelist = xrd_dirlist(skim_dir + skim_version)
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

def get_ntuple_filenames(
  dType, skim_version, years=years,
  batch=None, nbatches=None,
  redirector="ndcms.crc.nd.edu",
  test=False):
  '''Get input filenames'''
  #redirector = "deepthought.crc.nd.edu"
  #redirector = "cmsxrootd.fnal.gov"
  ntuple_location = f"root://{redirector}//store/group/lpcsusyphotons/{skim_version}"
  txtfile = f"{filelist_dir}/{skim_version}_filelist.txt"

  filenames = []
  with open(txtfile, "r") as f:
    for line in f:
      if tag_filename(dType, line, years=years):
        filenames.append(f"{ntuple_location}/{line[:-1]}")

  if test:
    # Write filenames to tmp file in cache
    tmp_file = f"{top_dir}/test/tmp/{dType}_{skim_version}_filelist.txt"
    with open(tmp_file, "w") as f:
      f.write("\n".join(filenames))


  if batch is not None and nbatches is not None:
    file_batches = np.array_split(filenames, nbatches)
    filenames = file_batches[batch]

  if not any(filenames):
    l.log(f"No files found for {dType} {skim_version} {years}", 0)
  
  return filenames


def characterize_skims(update_xrd=False, max_verbosity=0):
  # Get the number of files for each dType and year in each skim and store it in a pandas dataframe with columns: dType_year and rows: skim
  import pandas as pd
  
  tags = ["Skims", "TreeMaker","skims"]
  all_skims = xrd_dirlist(skim_dir) if update_xrd else [f.replace("_filelist.txt","") for f in os.listdir(filelist_dir)]
  skims = [ d for d in all_skims if any(tag in d for tag in tags) ]
  df = pd.DataFrame(index=skims, columns=['{}_{}'.format(dType, year) for dType in all_dTypes for year in years])
  for skim in skims:
    filelist = get_filelist(skim, max_verbosity=max_verbosity)
    for dType in all_dTypes:
      for year in years:
        df.loc[skim, '{}_{}'.format(dType, year)] = len([f for f in filelist if dType in f and year in f])

  # Save dataframe to YAML
  import yaml

  # Convert the dataframe to a dictionary
  data = df.transpose().to_dict()

  # Save the dictionary to a YAML file
  with open('skim_info.yaml', 'w') as f:
      yaml.dump(data, f, default_flow_style=False)

trigger_cache = f'{cache_dir}/trigger_index.json'
def get_trigger_index(dType, year, ntuple_version,
                      trigger_name='default', 
                      max_verbosity=1):
  
  l.set_max_verbosity(max_verbosity)

  if trigger_name == 'default':
    trigger_name = trigger_names[year]

  trigger_dict = {}
  cached = False
  if os.path.exists(trigger_cache):
    with open(trigger_cache, 'r') as f:
      trigger_dict = json.load(f)
      if (dType in trigger_dict.keys()) and \
          (ntuple_version in trigger_dict[dType].keys()) and \
          (year in trigger_dict[dType][ntuple_version].keys()):
        cached = True

  if not cached:
    l.log("Trigger index not cached, caching now...", 1)
    cache_trigger_indicies(ntuple_version, dTypes=[dType])
    with open(trigger_cache, 'r') as f:
      trigger_dict = json.load(f)

  trigger_index = trigger_dict[dType][ntuple_version][year]
  
  return trigger_index

def cache_trigger_indicies(ntuple_version,
                            dTypes=all_dTypes,
                           tree_name = 'TreeMaker2/PreSelection'):
  
  l.log('Caching trigger indicies for {} ntuple version {}'.format(dTypes, ntuple_version), 1)

  # Load cache if it exists
  if os.path.exists(trigger_cache):
    with open(trigger_cache, 'r') as f:
      trigger_dict = json.load(f)
  else:
    trigger_dict = {}

  for dType in dTypes:
    for year in years:
      filenames = get_ntuple_filenames(dType, ntuple_version, years=[year])
      if not any(filenames):
        l.log('No files found for {} in year {}'.format(dType, year), 2)
        continue
      l.log('Getting trigger index for {} in year {}'.format(dType, year), 2)
      filename = filenames[0]
      trigger_name = trigger_names[year]

      # Print trigger names
      output = get_tree_print(filename, tree_name, 0)
      
      # Parse the output
      parsed_names = re.findall(r'HLT_[\w\d_]*', output.decode('utf-8'))
      trigger_index = parsed_names.index(trigger_name)

      # Save to dictionary
      if dType not in trigger_dict.keys():
        trigger_dict[dType] = {}
      if ntuple_version not in trigger_dict[dType].keys():
        trigger_dict[dType][ntuple_version] = {}
      trigger_dict[dType][ntuple_version][year] = trigger_index

  # Write to cache
  with open(trigger_cache, 'w') as f:    
    json.dump(trigger_dict, f, indent=2)

def get_tree_print(filename, tree_name, batch):
    import subprocess
    # Write a temporary ROOT macro
    macro_dir = '/tmp/atownse2/'
    if not os.path.exists(macro_dir):
        os.makedirs(macro_dir)

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
  cache_trigger_indicies('skimsforAustin')
  cache_trigger_indicies('Skims_mediumPhotonV10')