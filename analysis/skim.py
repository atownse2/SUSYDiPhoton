import os
import sys

import concurrent.futures

from ROOT import TFile
from ROOT.TMath import Sqrt, Cos

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(top_dir)

from analysis.utils import sample_info as si
from analysis.utils import logger
from analysis.utils import condor


from analysis import selections
from analysis import outputs as out


l = logger.Logger()

all_input_versions = {
    "SR" : 'Skims_mediumPhotonV10',
    "EG" : 'Skims_mediumPhotonV10',
    "DY" : 'skimsforAustin',
}

signal_input_version = 'TreeMakerRandS_signal_fragmentedv9'

condor_shell_script = f"{top_dir}/analysis/utils/run_skim.sh"

job_splittings = {
    "data" : 200,
    "signal" : 500,
    "DYJetsToLL" : 20,
    "WGJets" : 20,
    "TTGJets" : 10,
}

def submit_skims(
        dType, analysis_region, git_tag,
        nbatches=None, 
        test=False,
        test_batch=False,
        verbosity=0):

    l.set_max_verbosity(verbosity)

    years = si.years
    if dType == "DYJetsToLL":
        years = ['2016','2017']
    elif dType == "signal":
        years = ['2016']

    if nbatches is None:
        nbatches = job_splittings[dType]
    
    if test or test_batch:
        nbatches = nbatches*20

    l.log(f"Submitting {nbatches} batches for {dType} in each {years}", 1)

    for year in years:
        for batch in range(nbatches):
            arguments = f"{dType} {analysis_region} {year} {nbatches} {batch} {verbosity}"
            
            tag = out.get_output_filename(dType, year, analysis_region, git_tag, "", "", tag_only=True)
            job_name = f"{tag}_{batch}"

            if test and not test_batch:
                condor.submit_local(condor_shell_script, arguments, job_name)
            else:
                condor.clear_logs(job_name)
                condor.submit_batch(condor_shell_script, arguments, job_name)
            
            if test or test_batch:
                return

class SkimTuples:

    def __init__(self,
                 dType, year, analysis_region, git_tag,
                 verbosity=0,
                 batch=None, 
                 nbatches=None):
        
        l.set_max_verbosity(verbosity)

        self.dType = dType
        self.analysis_region = analysis_region
        self.year = year
        self.git_tag = git_tag
        self.batch = batch
        self.nbatches = nbatches
        self.verbosity = verbosity

        self.isMC = dType != 'data'

    # When returning attributes look in task_config first
    def __getattr__(self, name):
        if hasattr(self.task_config, name):
            return getattr(self.task_config, name)
        else:
            return getattr(self, name)

    def run(self):
        self.outputs = out.Outputs(
            self.dType,
            self.year,
            self.analysis_region,
            self.git_tag,
            batch=self.batch)

        if self.dType == 'signal':
            self.input_version = signal_input_version
        else:
            self.input_version = all_input_versions[self.analysis_region]

        filenames = si.get_ntuple_filenames(
            self.dType,
            self.input_version,
            years=[self.year],
            batch=self.batch,
            nbatches=self.nbatches)
        
        l.log(f"Processing {len(filenames)} files", 1)
        self.process_files(filenames)

        self.outputs.save(verbosity=self.verbosity)

    def process_file(self, filename):
        entries = []
        if self.isMC:
            genelectrons = []

        l.log("Processing " + filename, 1)
        ### Setup
        dataFile = si.DataFile(filename)

        dataset = dataFile.dataset
        year = dataFile.year
        trigger_index = dataFile.trigger_index

        cutflow = self.outputs.cutflows[dataset]
        
        l.log(f"Dataset: {dataset} Year: {year} isMC: {self.isMC} Trigger Index: {trigger_index}", 2)

        ### Open file
        try:
            l.log("Opening File " + filename)
            f = TFile.Open(filename)
            t = f.Get('TreeMaker2/PreSelection')
        except:
            l.log(f"Could not open file {filename}", 0)
            return

        if self.isMC:
            counter = f.Get('tcounter')
            cutflow['total_mc'] += counter.GetEntries()

        l.log("Starting Event Loop",2)
        for event in t:

            if self.isMC:
                w = t.CrossSection*t.puWeight
            else:
                w = 1

            if self.isMC:
                # Gen electrons
                for iPar, genPar in enumerate(t.GenParticles):
                    if abs(t.GenParticles_PdgId[iPar]) != 11: continue
                    for iPho, pho in enumerate(t.Photons):
                        if selections.deltaR(pho, genPar) < 0.1:
                            genelectrons.append({
                                'genPt': genPar.Pt(),
                                'genEta': genPar.Eta(),
                                'genPhi': genPar.Phi(),
                                'hasPixelSeed': t.Photons_hasPixelSeed[iPho],
                                'met': t.HardMETPt,
                                'bdt': t.mva_BDT,
                                'nvtx': t.NVtx,
                                'weight': w,
                            })

            E = selections.EventSelection(t, trigger_index, self.analysis_region, cutflow)

            if not E.pass_selection: continue

            eg = E.candidates

            # If MinMt branch is not present, calculate it
            if hasattr(t, 'MinMt'):
                min_mT = t.MinMt
            else:
                # Calculate min mT
                mT = lambda pt, pt_miss, phi, phi_miss: Sqrt(2*pt*pt_miss*(1-Cos(phi-phi_miss)))
                if len(t.JetsAUX) == 0:
                    min_mT = None
                else:
                    min_mT = min([mT(j.Pt(), t.HardMETPt, j.Phi(), t.HardMETPhi) for j in t.JetsAUX])

            entry = {
                'dataset': dataset,
                'lead_pt': eg[0]['4vec'].Pt(),
                'subl_pt': eg[1]['4vec'].Pt(),
                'lead_eta': eg[0]['4vec'].Eta(),
                'subl_eta': eg[1]['4vec'].Eta(),
                'lead_photonWP': eg[0]['photonWP'],
                'subl_photonWP': eg[1]['photonWP'],
                'lead_hasPixelSeed': eg[0]['xseed'],
                'subl_hasPixelSeed': eg[1]['xseed'],
                "invariant_mass": (eg[0]['4vec'] + eg[1]['4vec']).M(),
                'min_mt': min_mT,
                'nvtx': t.NVtx,
                'met': t.HardMETPt,
                'bdt': t.mva_BDT,
                'weight': w,
                }
            
            if self.dType == 'signal':
                entry['SignalParameters'] =  [par for par in t.SignalParameters]

            # Can add gen info to entry here
            if self.isMC:
                selections.genMatching(t, eg)
                entry['lead_genMatch'] = eg[0]['genMatch']['name']
                entry['subl_genMatch'] = eg[1]['genMatch']['name']

                # Loop over gen electrons, match to photons, and record if they have a pixel seed or not
                n_gen_electrons = 0
                n_gen_electrons_matched_electron = 0
                n_gen_electrons_matched_photon = 0
                for iPar, genPar in enumerate(t.GenParticles):
                    if abs(t.GenParticles_PdgId[iPar]) != 11: continue
                    n_gen_electrons += 1
                    pho_match = None
                    for iPho, pho in enumerate(t.Photons):
                        if selections.deltaR(pho, genPar) < 0.1:
                            pho_match = iPho
                            break
                    
                    if pho_match is not None:
                        if t.Photons_hasPixelSeed[pho_match]:
                            n_gen_electrons_matched_electron += 1
                        else:
                            n_gen_electrons_matched_photon += 1
                
                entry['n_gen_electrons'] = n_gen_electrons
                entry['n_gen_electrons_matched_electron'] = n_gen_electrons_matched_electron
                entry['n_gen_electrons_matched_photon'] = n_gen_electrons_matched_photon

            entries.append(entry)

        # After event loop
        f.Close()

        outputs = {'events': entries}
        if self.isMC:
            outputs['genelectrons'] = genelectrons

        return outputs

    def process_files(self, filenames):
        # Run multithreaded

        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = executor.map(self.process_file, filenames)
        
        for result in results:
            self.outputs.events.fill(result['events'])
            if self.isMC:
                self.outputs.genelectrons.fill(result['genelectrons'])

if __name__ == "__main__":
    import time

    from analysis.utils import version
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dType", help="Data type to process")
    parser.add_argument("analysis_region", help="Analysis region to process")
    parser.add_argument("--submit_batches", action="store_true", help="Submit batch jobs")
    parser.add_argument("--git_tag", "-tag", type=str, default=version.get_last_tag(), help="Git tag to use for output")
    parser.add_argument("-y", "--year", type=str, default=None, help='years to run')
    parser.add_argument("-n", "--nbatches", type=int, default=1, help="Number of batches")
    parser.add_argument("-b", "--batch", type=int, default=0, help="Batch number")
    parser.add_argument("-v", "--verbosity", type=int, default=0, help="Verbosity level")
    parser.add_argument("-t", "--test", action="store_true", help="Run test")
    parser.add_argument("--test_batch", action="store_true", help="Run test batch")
    args = parser.parse_args()

    if args.test:
        args.year = "2016"
        args.nbatches = 10000
        args.verbosity = 2
        out.output_dir = f"{top_dir}/test/outputs"
        out.batch_output_dir = f"{top_dir}/test/outputs/batch"

    if args.submit_batches or args.test_batch:
        version.check_git()
        args.git_tag = version.get_last_tag()
        submit_skims(args.dType, args.analysis_region, args.git_tag, test=args.test, verbosity=args.verbosity, test_batch=args.test_batch)
    else:
        if args.year is None:
            raise ValueError("Must specify year when not submitting batches")
        
        t1 = time.time()

        SkimTuples(
            args.dType,
            args.year,
            args.analysis_region,
            args.git_tag,
            verbosity=args.verbosity,
            batch=args.batch, nbatches=args.nbatches).run()
        t2 = time.time()
        print(f"Time elapsed: {t2-t1:.2f} s")


