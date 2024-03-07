import os
import sys

import json

import numpy as np

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(top_dir)

import analysis.fitting as fit
from analysis import outputs as out


tf_cache = f"{top_dir}/cache/tf"
if not os.path.exists(tf_cache):
    os.makedirs(tf_cache)


def compare_tfs(tf1, tf2):
    assert tf1.years == tf2.years
    for year in tf1.years:
        print(f"{year}:")
        print(f"{tf1.name}: tf : {tf1.tf[year]}, tf err: {tf1.tf_err[year]}")
        print(f"{tf2.name}: tf : {tf2.tf[year]}, tf err: {tf2.tf_err[year]}")

def round_tf(tf, tf_err, n=3):
    # Find the number of decimal places to the first two significant figures
    # ndec = 0
    # for i, s in enumerate(str(tf_err).split('.')[1]):
    #     if s != '0':
    #         ndec = i
    #         break
    
    # tf = round(tf, ndec+n)
    # tf_err = round(tf_err, ndec+n)

    return tf, tf_err


class BaseTF:
    def __init__(self, isMC, git_tag, years=None, remake=False):

        self.isMC = isMC
        self.git_tag = git_tag

        self.dType = "DYJetsToLL" if isMC else "data"

        self.years = ["2016", "2017", "2018"] if years is None else years
        if isMC:
            self.years.remove("2018")

        self.info = {year: None for year in self.years}
        self.cache = f"{tf_cache}/{self.name}_{self.dType}.json"

        self.init_tf(remake)

    def init_tf(self, remake):
        if os.path.exists(self.cache) and not remake:
            self.load_tf()
        else:
            self.calculate_tf()
            self.save_tf()

    @property
    def tf(self):
        if self.isMC:
            self.info["2018"] = self.info["2017"]
        return {year: info["tf"] for year, info in self.info.items()}
    
    @property
    def tf_err(self):
        if self.isMC:
            self.info["2018"] = self.info["2017"]
        return {year: info["tf_err"] for year, info in self.info.items()}

    def load_tf(self):
        with open(self.cache, 'r') as f:
            self.info = json.load(f)
    
    def save_tf(self):
        with open(self.cache, 'w') as f:
            json.dump(self.info, f, indent=4)

class CombinedTF(BaseTF):

    def __init__(self, isMC, git_tag, **kwargs):
        self.name = "CombinedTF"
        super().__init__(isMC, git_tag, **kwargs)

    def calculate_tf(self):
        for year in self.years:
            self.info[year] = {"tf": {}, "tf_err": {}}

            df = out.load_outputs(self.dType, year, "DY", self.git_tag).events.df

            BB = (abs(df['lead_eta'])<1.4442) & (abs(df['subl_eta'])<1.4442)
            BE = (abs(df['lead_eta'])<1.4442) & (abs(df['subl_eta'])>1.566)
            EB = (abs(df['lead_eta'])>1.566) & (abs(df['subl_eta'])<1.4442)
            EE = (abs(df['lead_eta'])>1.566) & (abs(df['subl_eta'])>1.566)

            ee = df['lead_hasPixelSeed'] & df['subl_hasPixelSeed']
            eg = df['lead_hasPixelSeed'] & ~df['subl_hasPixelSeed']
            ge = ~df['lead_hasPixelSeed'] & df['subl_hasPixelSeed']

            # Scale the eBeB and eEeE weights by 2
            df.loc[ee&BB,'weight']*=2
            df.loc[ee&EE,'weight']*=2

            regions = {
                'barrel' : {
                    'numerator' : (BB & (eg | ge)) | (BE & ge) | (EB & eg),
                    'denominator' : (ee & BB) | (ee & BE) | (ee & EB)},
                'endcap' : {
                    'numerator' : (EE & (eg | ge)) | (BE & eg) | (EB & ge),
                    'denominator' : (ee & EE) | (ee & BE) | (ee & EB)}
            }

            for region, cuts in regions.items():
                num_fit = fit.fit(df[cuts['numerator']], f'{self.dType} {region} numerator {year}', 'binned')
                den_fit = fit.fit(df[cuts['denominator']], f'{self.dType} {region} denominator {year}', 'binned')

                tf = num_fit['nsig']/den_fit['nsig']
                tf_err = np.sqrt( (num_fit['nsig_err']/den_fit['nsig'])**2 + (den_fit['nsig_err']*num_fit['nsig']/(den_fit['nsig']**2))**2 )

                tf, tf_err = round_tf(tf, tf_err)

                self.info[year]["tf"][region] = tf
                self.info[year]["tf_err"][region] = tf_err

class NaiveTF(BaseTF):
    
        def __init__(self, isMC, git_tag, **kwargs):
            self.name = "NaiveTF"
            super().__init__(isMC, git_tag, **kwargs)
    
        def calculate_tf(self):
            for year in self.years:            
                df = out.load_outputs(self.dType, year, "DY", self.git_tag).events.df
    
                ee_CR = df['lead_hasPixelSeed'] & df['subl_hasPixelSeed']
                eg_CR = df['lead_hasPixelSeed'] != df['subl_hasPixelSeed']
    
                # Provide estimates for the transfer factor using
                same_regions = {
                    "BB" : (abs(df['lead_eta'])<1.4442) & (abs(df['subl_eta'])<1.4442),
                    "EE" : (abs(df['lead_eta'])>1.566) & (abs(df['subl_eta'])>1.566)
                }

                self.info[year] = {"tf": {}, "tf_err": {}}
                for region, cut in same_regions.items():
                    ee_fit = fit.fit(df[ee_CR & cut], f'{self.dType} ee CR {year} {region}', 'binned')
                    eg_fit = fit.fit(df[eg_CR & cut], f'{self.dType} eg CR {year} {region}', 'binned')
    
                    tf = eg_fit['nsig']/(2*ee_fit['nsig'])
                    tf_err = np.sqrt( (eg_fit['nsig_err']/(2*ee_fit['nsig']))**2 + (ee_fit['nsig_err']*eg_fit['nsig']/(2*ee_fit['nsig']**2))**2 )
    
                    tf, tf_err = round_tf(tf, tf_err)

                    self.info[year]["tf"][region] = tf
                    self.info[year]["tf_err"][region] = tf_err

class SimpleTF(BaseTF):

    def __init__(self, isMC, git_tag, **kwargs):
        self.name = "SimpleTF"
        super().__init__(isMC, git_tag, **kwargs)

    def calculate_tf(self):
        for year in self.years:            
            df = out.load_outputs(self.dType, year, "DY")[0]

            ee_CR = df['lead_hasPixelSeed'] & df['subl_hasPixelSeed']
            eg_CR = df['lead_hasPixelSeed'] != df['subl_hasPixelSeed']

            ee_fit = fit.fit(df[ee_CR], f'{self.dType} ee CR {year}', 'binned')
            eg_fit = fit.fit(df[eg_CR], f'{self.dType} eg CR {year}', 'binned')

            tf = eg_fit['nsig']/(2*ee_fit['nsig'])
            tf_err = np.sqrt( (eg_fit['nsig_err']/(2*ee_fit['nsig']))**2 + (ee_fit['nsig_err']*eg_fit['nsig']/(2*ee_fit['nsig']**2))**2 )

            tf, tf_err = round_tf(tf, tf_err)

            self.info[year] = {"tf": tf, "tf_err": tf_err}