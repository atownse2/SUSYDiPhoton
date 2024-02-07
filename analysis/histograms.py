from array import array
import numpy as np

import hist

from analysis.utils import sample_info as si

final_states = [
    'eBeB', 'eBeE', 'eEeB', 'eEeE', # ee
    'gBgB', 'gBgE', 'gEgB', 'gEgE', # gg
    'eBgB', 'gBeB', #eg BB
    'eEgE', 'gEeE', #eg EE
    'eBgE', 'gEeB', 'eEgB', 'gBeE'] # eg EB

reco_states = ['gg', 'eg', 'ge', 'ee']

gg_final_states = [fs for fs in final_states if fs.count('g') == 2]
eg_final_states = [fs for fs in final_states if fs.count('g') == 1 and fs.count('e') == 1]

det_regions = ['barrel', 'endcap']

final_state_tf_region_map = {
    'eBgB' : 'barrel',
    'gBeB' : 'barrel',
    'eBgE' : 'barrel',
    'gBeE' : 'endcap',
    'eEgB' : 'endcap',
    'gEeB' : 'barrel'
}

prediction_map = { 'gBgB' : ['eBgB', 'gBeB'],
                   'gBgE' : ['eBgE', 'gBeE'], 
                   'gEgB' : ['eEgB', 'gEeB']}

gg_regions = ['BB', 'BE', 'EB']

###

genEvents = ['genGJet', 'genEG']
genTypes = ["genEle", "genMu", "genTau", "genPho", "genJet"]
gg_cutflow = ['all', 'fakehardMET', 'realhardMET', 'realhardMET+lostLepton']

egtf_bins = {'pt' : array("d", [80,100,150,200]),
             'region': array("d", [0,1]), # Barrel = 1, Endcap = 0
             'eta' : array("d", [0,0.4, 0.8, 1.2, 1.4442, 1.566, 2.0, 2.5]),
             'nvtx' : array("d", np.arange(0,64,8)),
             'invariant_mass' : array("d", np.linspace(60,120, num=32)),
             }

analysis_bins = {'bdt' : array("d",[-1, -0.13, 0.03, 1]),
                 'met' : array("d",[130, 150, 185, 250]),
}

var_bins = {'bdt' : array("d", np.linspace(-1,1, num=32)),
        'pt' : array("d", [80,100,150,200,300]),
        'eta': array("d", [0, 0.8, 1.4442, 1.566, 2.0, 2.5]),
        'met' : array("d", [130, 150, 185, 250]),
        'invariant_mass' : array("d", np.linspace(60,120, num=32)),
        'nvtx' : array("d", range(0, 64)),
}

var_labels = { 'lead_pt' : 'Leading (Fake) Photon $p_{T}$ [GeV]',
              'subl_pt' : 'Sub-Leading (Fake) Photon $p_{T}$ [GeV]',
              'lead_eta' : 'Leading (Fake) Photon $\eta$',
              'subl_eta' : 'Sub-Leading (Fake) Photon $\eta$',
              'nvtx' : '$N_{vtx}$',
              'invariant_mass' : 'Invariant Mass [GeV]',
              'bdt' : 'BDT Score',
              'met' : 'hardMET [GeV]',}


def get_var_bin_names(
    var,
    underflow=None, overflow=None,
    use_egtf_bins=False,
    use_analysis_bins = False,
    use_bdt_bins=False):
    if use_bdt_bins:
        return ['lowBDT', 'medBDT', 'highBDT']
    elif use_egtf_bins:
        bins = egtf_bins[var]
    elif use_analysis_bins:
        bins = analysis_bins[var]
    else:
        bins = var_bins[var]
    var_bin_names = [f"{var}{bins[i]}to{bins[i+1]}" for i in range(len(bins)-1)]
    
    if underflow is not None:
        var_bin_names = [f"{var}{underflow}to{bins[0]}"] + var_bin_names
    if overflow is not None:
        var_bin_names = var_bin_names + [f"{var}{bins[-1]}to{overflow}"]

    return var_bin_names

def get_bdt_bin_name(bdtVal):
    bin_names = get_var_bin_names('bdt', use_bdt_bins=True)
    for i in range(len(analysis_bins['bdt'])-1):
        if bdtVal <= analysis_bins['bdt'][i+1]:
          return bin_names[i]
    return None

# Hist axes
dataset_axis = hist.axis.StrCategory( [], name="dataset", label="Dataset", growth=True)
year_axis = hist.axis.StrCategory( si.years, name="year", label="Year")
lead_gen_match_axis = hist.axis.StrCategory([], name="lead_gen_match", label="Lead Gen Match", growth=True)
subl_gen_match_axis = hist.axis.StrCategory([], name="subl_gen_match", label="Trail Gen Match", growth=True)
lead_pt_axis = hist.axis.Variable(var_bins['pt'], name="lead_pt", label="Lead Photon p_{T} [GeV]")
subl_pt_axis = hist.axis.Variable(var_bins['pt'], name="subl_pt", label="Trail Photon p_{T} [GeV]")
lead_eta_axis = hist.axis.Variable(var_bins['eta'], name="lead_eta", label="Lead Photon #eta")
subl_eta_axis = hist.axis.Variable(var_bins['eta'], name="subl_eta", label="Trail Photon #eta")
invariant_mass_axis = hist.axis.Variable(var_bins['invariant_mass'], name="invariant_mass", label="Invariant Mass [GeV]")

ar_axis = hist.axis.StrCategory([], name="analysis_region", label="Analysis Region", growth=True)
final_state_axis = hist.axis.StrCategory([], name="final_state", label="Final State", growth=True) 
met_axis = hist.axis.Variable(var_bins['met'], name="met", label="MET [GeV]")
bdt_axis = hist.axis.Variable(var_bins['bdt'], name="bdt", label="BDT Score")
nvtx_axis = hist.axis.Variable(egtf_bins['nvtx'], name="nvtx", label="N_{vtx}", underflow=False)

egtf_axes = {'lead_pt': hist.axis.Variable(egtf_bins['pt'], name="lead_pt", label="Leading Photon p_{T} [GeV]"),
             'subl_pt': hist.axis.Variable(egtf_bins['pt'], name="subl_pt", label="Sub-Leading Photon p_{T} [GeV]"),
             'lead_eta': hist.axis.Variable(egtf_bins['eta'], name="lead_eta", label="Leading Photon \eta"),
             'subl_eta': hist.axis.Variable(egtf_bins['eta'], name="subl_eta", label="Sub-Leading Photon \eta"),
             'nvtx': hist.axis.Variable(egtf_bins['nvtx'], name="nvtx", label="N_{vtx}", underflow=False),
             'invariant_mass': hist.axis.Variable(egtf_bins['invariant_mass'], name="invariant_mass", label="Invariant Mass [GeV]"),
}
analysis_axes = {'bdt': hist.axis.Variable(analysis_bins['bdt'], name="bdt", label="BDT Score", underflow=False, overflow=False),
                 'met': hist.axis.Variable(analysis_bins['met'], name="met", label="MET [GeV]", underflow=False, overflow=True),}

class SimpleHistogram:
    def __init__(self,):
        self.values = None
        self.variances = None
        self.bins = None

    def from_numpy(self, np_hist, variance=None):
        self.values = np_hist[0]
        self.bins = np_hist[1:]
        self.variances = variance
        return self
    def from_hist(self, hist, flow=True):
        np_hist = hist.to_numpy(flow=flow)
        self.values = np_hist[0]
        self.bins = np_hist[1:]
        self.variances = hist.variances(flow=flow)
        return self
    def from_self(self, values, variances, bins):
        self.values = values
        self.variances = variances
        self.bins = bins
        return self

    def to_numpy(self):
        return (self.values, *self.bins)

    def same_bins(self, other):
        return all([np.array_equal(self.bins[i], other.bins[i]) for i in range(len(self.bins))])

    def __add__(self, other):
        if self.values is None:
            return other

        assert self.same_bins(other)
        
        self.values += other.values
        self.variances += other.variances
        return self

    def __truediv__(self, other):
        assert self.same_bins(other)
        assert np.all(other.values != 0)

        values = self.values/other.values
        variances = (self.variances/other.values**2 + self.values**2*other.variances/other.values**4)

        return SimpleHistogram().from_self(values, variances, self.bins)

    def scale(self, factor, err):
        self.values *= factor
        self.variances = (self.values*err)**2 + self.variances*factor**2
        return self

def analysis_hist(data):
    return hist.Hist(
        analysis_axes['bdt'],
        analysis_axes['met'],
        storage=hist.storage.Weight(),
        ).fill(
        bdt=data['bdt'],
        met=data['met'],
        weight=data['weight'])

def invariantMass(data, hist_name):
    return hist.Hist(
        invariant_mass_axis,
        storage=hist.storage.Weight()
        ).fill(
        invariant_mass=data['invariantMass'],
        weight=data['weight'],
        )

def genBackground(data):
    genMatches = ["genPho", "genEle", "genMu", "genTau", "genJet"]
    return hist.Hist( 
        hist.axis.StrCategory(genMatches, name="lead_genMatch", label="Lead Photon Gen Match"),
        hist.axis.StrCategory(genMatches, name="subl_genMatch", label="Sublead Photon Gen Match"),
        storage=hist.storage.Weight()
        ).fill(
        lead_genMatch=data['lead_genMatch'],
        subl_genMatch=data['subl_genMatch'],
        weight=data['weight'],
        )

def genFake():
    return hist.Hist(
        dataset_axis,
        hist.axis.Regular(50, 0, 200, name="pt", label="Gen Electron p_{T} [GeV]"),
        hist.axis.Regular([0, 0.4, 0.8, 1.4442, 1.566, 2.0, 2.5], name="eta", label="Gen Electron #eta"),
        hist.axis.Boolean(name="hasPixelSeed", label="Gen Electron hasPixelSeed"),
        bdt_axis,
        met_axis,

        storage=hist.storage.Weight()
        )