import os

import git

import pandas as pd

from analysis.utils import sample_info as si
from analysis.utils import logger
from analysis.utils import condor
from analysis.utils import version

from analysis import skim

l = logger.Logger()

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

output_dir = f"{top_dir}/outputs"
batch_output_dir = f"{top_dir}/outputs/batch"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(batch_output_dir):
    Warning(f"Creating batch output directory {batch_output_dir} for the first time, please set the correct permissions for condor.")
    os.makedirs(batch_output_dir)


def get_output_filename(dType, year, analysis_region, git_tag, extension, output_type,
                       batch=None, tag_only=False):

    tag = f'{dType}_{year}_{analysis_region}_{git_tag}'
    if tag_only:
        return tag
    
    if batch is None:
        outdir = f"{output_dir}/{output_type}"
    else:
        outdir = f"{batch_output_dir}/{output_type}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if batch is not None:
        tag = f"{tag}_{batch}"

    filepath = f"{outdir}/{tag}.{extension}"

    return filepath

def load_outputs(dType, era, analysis_region, git_tag):
    return Outputs(dType, era, analysis_region, git_tag).load()

class Outputs:
    def __init__(self, dType, era, analysis_region, git_tag, batch=None, test_load=False):

        self.dType = dType
        self.era = era
        self.years = si.get_years_from_era(era)
        self.analysis_region = analysis_region
        self.git_tag = git_tag

        self.test_load = test_load
        self.load_instance = False

        self.types = [
            Cutflows(dType, era, analysis_region, git_tag, batch=batch),
            Events(dType, era, analysis_region, git_tag, batch=batch),]

        if dType != 'data':
            self.types.append(GenElectrons(dType, era, analysis_region, git_tag, batch=batch))
        
        self.outputs = { output.output_type: output for output in self.types}
    
    def __getattr__(self, name):
        if name in self.outputs:
            if not self.outputs[name].loaded and self.load_instance:
                self.load(lazy=False, test=self.test_load, output_names=name)
            return self.outputs[name]

    def __add__(self, other):
        for output in self.outputs.values():
            output += other.outputs[output.output_type]
        return self

    def load(self, lazy=True, output_names=None, test=False):
        if lazy:
            self.load_instance = True
            return self

        if output_names is None:
            output_names = self.outputs.keys()
        if isinstance(output_names, str):
            output_names = [output_names]
        
        for output_name in output_names:
            self.outputs[output_name].scaled = True
            output_type = type(self.outputs[output_name])
            for year in self.years:
                self.outputs[output_name] += output_type(self.dType,
                                                         year,
                                                         self.analysis_region,
                                                         self.git_tag).load(test=test)
            self.outputs[output_name].loaded = True
        
        return self
    
    def save(self):
        for output in self.outputs.values():
            output.save()

class Output:
    def __init__(self, dType, year, analysis_region, git_tag, batch=None):
        self.dType = dType
        self.year = year
        self.analysis_region = analysis_region
        self.git_tag = git_tag
        self.batch = batch

        self.loaded = False
        self.scaled = False

    @property
    def filename(self):
        return self.get_filename()

    def get_filename(self, tag_only=False, batch=None):
        if batch is None: batch = self.batch
        return get_output_filename(
            self.dType, self.year, self.analysis_region, self.git_tag, self.extension, self.output_type, batch=batch, tag_only=tag_only)
    
    def can_add(self, other):
        if self.dType != other.dType: return False
        if self.analysis_region != other.analysis_region: return False
        if not (self.scaled and other.scaled) and self.year != other.year: return False
        return True

    def load(self, test=False):
        if not os.path.exists(self.filename):
            # Look for batch files
            batch_dir = os.path.dirname(self.get_filename(batch=0))
            file_tag = self.get_filename(tag_only=True)
            print(f"Looking for batch files in {batch_dir} with tag {file_tag}")
            batch_filenames = [f"{batch_dir}/{f}" for f in os.listdir(batch_dir) if f.startswith(file_tag) and f.endswith(self.extension)]

            if len(batch_filenames) == 0:
                print(f"Batch events for {file_tag} do not exist")
                return None

            # Check the condor logs for errors
            condor.check_logs(file_tag)

            # Check for missing batch files
            if len(batch_filenames) != skim.job_splittings[self.dType]:
                print(f"Warning: {len(batch_filenames)} batch files found for {file_tag}, expected {skim.job_splittings[self.dType]}")
                input("Press Enter to continue...")
            
            for batch_filename in batch_filenames:
                self += type(self)(self.dType,
                                   self.year,
                                   self.analysis_region,
                                   self.git_tag).load_from_file(batch_filename)

            # Remove batch files
            if not test:
                for batch_filename in batch_filenames:
                    os.remove(batch_filename)
            
        else:
            self.load_from_file(self.filename)
        
        self.loaded = True
        self.save() # Save combined file

        if self.dType != 'data' and hasattr(self, 'scale_mc'):
            self.scale_mc()

        return self

class Cutflow():
    def __init__(self, dType, year, analysis_region):
        self.dType = dType
        self.year = year
        self.analysis_region = analysis_region

        self.cutflow = {}
    
    def __add__(self, other):
        for cut, count in other.cutflow.items():
            if cut not in self.cutflow:
                self.cutflow[cut] = count
            else:
                self.cutflow[cut] += count
        return self
    
    def keys(self):
        return self.cutflow.keys()

    def values(self):
        return self.cutflow.values()

    def items(self):
        return self.cutflow.items()

    def __getitem__(self, cut):
        if cut not in self.cutflow:
            self.cutflow[cut] = 0
        return self.cutflow[cut]

    def __setitem__(self, cut, count):
        self.cutflow[cut] = count
    
    def to_dict(self):
        return self.cutflow
    def from_dict(self, cutflow):
        self.cutflow = cutflow
        return self

class Cutflows(Output):
    def __init__(self, dType, year, analysis_region, git_tag, batch=None):
        super().__init__(dType, year, analysis_region, git_tag, batch=batch)

        self.info = {}

        self.output_type = 'cutflows'
        self.extension = 'json'
    
    @property
    def datasets(self):
        return self.info.keys()

    def print(self):
        for dataset, cutflow in self.info.items():
            print(f"{dataset}:")
            for cut, count in cutflow.cutflow.items():
                print(f"  {cut}: {count}")

    def __add__(self, other):
        assert self.can_add(other), "Cannot add Cutflows objects with different datasets"
        for dataset, cutflow in other.info.items():
            if dataset not in self.info:
                self.info[dataset] = cutflow
            else:
                self.info[dataset] += cutflow
        return self
    
    def __getitem__(self, dataset):
        if dataset not in self.info:
            self.info[dataset] = Cutflow(self.dType, self.year, self.analysis_region)
        return self.info[dataset]
    
    def save(self):
        import json
        to_save = {dataset: cutflow.to_dict() for dataset, cutflow in self.info.items()}
        with open(self.filename, 'w') as f:
            json.dump(to_save, f, indent=4)
    
    def load_from_file(self, filename):
        import json
        with open(filename, 'r') as f:
            from_load = json.load(f)
        for dataset, cutflow in from_load.items():
            self.info[dataset] = Cutflow(self.dType, self.year, self.analysis_region).from_dict(cutflow)

        return self

    def scale_mc(self):
        cutflow = Cutflow(self.dType, self.year, self.analysis_region)
        for dataset, info in self.info.items():
            total_mc = info['total_mc']
            for cut, counts in info.items():
                if cut == 'total_mc': continue
                cutflow[cut] += counts*1000*si.lumis[self.year]/total_mc
            cutflow['total_mc'] += total_mc

        self.scaled = True
        self.info['weighted_sum'] = cutflow
        return self

class Events(Output):

    def __init__(self, dType, year, analysis_region, git_tag, batch=None):
        super().__init__(dType, year, analysis_region, git_tag, batch=batch)

        self.df = None

        self.output_type = 'events'
        self.extension = 'parquet'
    
    def __add__(self, other):
        assert self.can_add(other), "Cannot add Events objects with different datasets"
        if self.df is None:
            self.df = other.df
        else:
            self.df = pd.concat([self.df, other.df])
        return self

    def fill(self, entries):
        if self.df is None:
            self.df = pd.DataFrame(entries)
        else:
            self.df = pd.concat([self.df, pd.DataFrame(entries)])
        
        return self

    def save(self):
        self.df.to_parquet(self.filename)
    
    def load_from_file(self, filename):
        self.df = pd.read_parquet(filename)
        return self

    def scale_mc(self):
        if self.dType == 'data':
            raise ValueError("Cannot scale data")

        # Get dataset info
        cutflow = Cutflows(self.dType, self.year, self.analysis_region, self.git_tag, batch=self.batch).load()

        for dataset in cutflow.datasets:
            total_mc = cutflow[dataset]['total_mc']
            if total_mc == 0:
                raise ValueError(f"Dataset {dataset} has zero events")
            self.df.loc[self.df.dataset == dataset, 'weight'] *= 1000*si.lumis[self.year]/total_mc
        
        # Remove dataset column
        self.df.drop(columns=['dataset'], inplace=True)

        self.scaled = True
        return self

class GenElectrons(Events):
    def __init__(self, dType, year, analysis_region, git_tag, batch=None):
        super().__init__(dType, year, analysis_region, git_tag, batch=batch)
        self.output_type = 'genelectrons'

class Histograms(Output):

    def __init__(self, dType, year, analysis_region, git_tag, batch=None):
        super().__init__(dType, year, analysis_region, git_tag, batch=batch)

        self.histograms = {}

        self.output_type = 'histograms'
        self.extension = 'pkl'
    
    def __add__(self, other):
        assert self.can_add(other), "Cannot add Histograms objects with different datasets"
        for hist in other.histograms:
            if hist not in self.histograms:
                self.histograms[hist] = other.histograms[hist]
            else:
                self.histograms[hist] += other.histograms[hist]

    def __getitem__(self, histname):
        return self.histograms[histname]

    def __setitem__(self, histname, hist):
        self.histograms[histname] = hist

    def save(self):
        import pickle
        with open(self.filename, 'wb') as f:
            pickle.dump(self.histograms, f)

    def load_from_file(self, filename):
        import pickle
        with open(filename, 'rb') as f:
            self.histograms = pickle.load(f)
        return self

    def scale_mc(self):
        if self.dType == 'data':
            raise ValueError("Cannot scale data")

        # Get dataset info
        cutflow = Cutflows(self.dType, self.year, self.analysis_region, self.git_tag, batch=self.batch).load()

        for histname, hist in self.histograms.items():
            for dataset in cutflow.datasets:
                total_mc = cutflow[dataset]['total_mc']
                if total_mc == 0:
                    raise ValueError(f"Dataset {dataset} has zero events")
                
                hist[{'dataset':dataset}]*=1000*si.lumis[self.year]/total_mc
        
            # Sum over datasets
            self.histograms[histname] = hist.sum('dataset')
        
        self.scaled = True
        return self
