import os
import sys

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(top_dir)

from analysis.utils import sample_info as si
from analysis.utils import logger
from analysis.utils import condor

from analysis import selections
from analysis import outputs as out

l = logger.Logger()

tuple_version = {
    "SR" : 'Skims_mediumPhotonV10',
    "EG" : 'Skims_mediumPhotonV10',
    "DY" : 'skimsforAustin',
}

condor_shell_script = f"{top_dir}/analysis/utils/run_skim.sh"

job_splittings = {
    'data' : 100,
    'WGJets': 15,
    'DYJetsToLL' : 20,
    'TTGJets' : 10,
}

def print_diff(repo):
    import subprocess
    diff_index = repo.index.diff(None)
    for diff in diff_index.iter_change_type('M'):  # 'M' for modified files
        if diff.a_path.endswith('.ipynb'): continue
        print(subprocess.check_output(['git', 'diff', '--', diff.a_path]).decode())

def check_git():
    import git
    repo = git.Repo(top_dir)

    if not repo.is_dirty(): return

    print("Repository has uncommitted changes")
    print_diff(repo)
    _ = input("Would you like to commit these changes now? y/n: ")
    if _ in ['n', 'N']: return

    repo.git.add(update=True)
    message = input("Enter commit message: ")
    repo.index.commit(message)


def submit_skims(dType, analysis_region,
        nbatches=None, 
        test=False,
        test_batch=False,
        verbosity=0):

    l.set_max_verbosity(verbosity)

    years = si.years
    if dType == "DYJetsToLL":
        years = ['2016','2017']

    if nbatches is None:
        nbatches = job_splittings[dType]
    
    if test or test_batch:
        nbatches = 1000

    l.log(f"Submitting {nbatches} batches for {dType} in each {years}", 1)

    for year in years:
        for batch in range(nbatches):
            arguments = f"{dType} {analysis_region} {year} {nbatches} {batch} {verbosity}"
            
            tag = out.get_output_filename(dType, year, analysis_region, "", "", tag_only=True)
            job_name = f"{tag}_{batch}"

            if test and not test_batch:
                condor.submit_local(condor_shell_script, arguments, job_name)
            else:
                condor.clear_logs(job_name)
                condor.submit_batch(condor_shell_script, arguments, job_name)
            
            if test or test_batch:
                return

class SkimTuples:

    def __init__(self, dType, analysis_region, year,
                 verbosity=0,
                 batch=None, 
                 nbatches=None):
        
        l.set_max_verbosity(verbosity)

        self.dType = dType
        self.isMC = dType != 'data'
        self.analysis_region = analysis_region
        self.year = year
        self.batch = batch
        self.nbatches = nbatches

    def run(self):
        self.outputs = out.Outputs(self.dType, self.year, self.analysis_region, batch=self.batch)

        filenames = si.get_ntuple_filenames(
            self.dType, tuple_version[self.analysis_region], years=[self.year],
            batch=self.batch, nbatches=self.nbatches)
        
        l.log(f"Processing {len(filenames)} files", 1)
        self.fill_outputs(filenames)

        self.outputs.save()

    def fill_outputs(self, filenames):
        '''Loop over files in filenames, fill outputs'''
        from ROOT import TFile
        from ROOT.TMath import Sqrt, Cos

        entries = []
        if self.isMC:
            genelectrons = []

        for filename in filenames:
            l.log("Processing " + filename, 1)
            ### Setup
            dataFile = si.DataFile(filename)

            dataset = dataFile.dataset
            year = dataFile.year
            trigger_index = dataFile.trigger_index

            cutflow = self.outputs.cutflows[dataset]
            
            l.log(f"Dataset: {dataset} Year: {year} isMC: {self.isMC} Trigger Index: {trigger_index}", 2)

            ### Open file
            l.log("Opening File " + filename)
            f = TFile.Open(filename)
            t = f.Get('TreeMaker2/PreSelection')

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
                    'lead_hasPixelSeed': eg[0]['xseed'],
                    'subl_hasPixelSeed': eg[1]['xseed'],
                    "invariant_mass": (eg[0]['4vec'] + eg[1]['4vec']).M(),
                    'min_mt': min_mT,
                    'nvtx': t.NVtx,
                    'met': t.HardMETPt,
                    'bdt': t.mva_BDT,
                    'weight': w,
                    }
                
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
        self.outputs.events.fill(entries)

        if self.isMC:
            self.outputs.genelectrons.fill(genelectrons)

if __name__ == "__main__":
    import time

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dType", help="Data type to process")
    parser.add_argument("analysis_region", help="Analysis region to process")
    parser.add_argument("--submit_batches", action="store_true", help="Submit batch jobs")
    parser.add_argument("-y", "--year", type=str, default=None, help='years to run')
    parser.add_argument("-n", "--nbatches", type=int, default=1, help="Number of batches")
    parser.add_argument("-b", "--batch", type=int, default=0, help="Batch number")
    parser.add_argument("-v", "--verbosity", type=int, default=0, help="Verbosity level")
    parser.add_argument("-t", "--test", action="store_true", help="Run test")
    parser.add_argument("--test_batch", action="store_true", help="Run test batch")
    args = parser.parse_args()

    if args.test:
        args.year = "2016"
        args.nbatches = 1000
        args.verbosity = 2
        out.output_dir = f"{top_dir}/test/outputs"
        out.batch_output_dir = f"{top_dir}/test/outputs/batch"

    if args.submit_batches or args.test_batch:
        check_git()
        submit_skims(args.dType, args.analysis_region, test=args.test, verbosity=args.verbosity, test_batch=args.test_batch)
    else:
        if args.year is None:
            raise ValueError("Must specify year when not submitting batches")
        
        t1 = time.time()
        SkimTuples(args.dType, args.analysis_region,
                  year=args.year,
                  verbosity=args.verbosity,
                  batch=args.batch, nbatches=args.nbatches).run()
        t2 = time.time()
        print(f"Time elapsed: {t2-t1:.2f} s")
