import sys

top = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton"
sys.path.append(top)
from analysis.metadata.sample_info import *
from analysis.hists import *
from analysis.selections import *

from analysis.logger import Logger

l = Logger()

### Configuration
photonWP = 'loose'
version = "101023v1"
###
ntuple_version  = 'skimsforAustin'
skimmer = 'egtf'
###

def make_hists(dType, subsets=True):
  ### Initialize histograms
  hists = OrderedDict()
  #For E->G fakerate prediction
  isMC = dType != 'data'
  
  for era in eras:
    for fs in final_states:
      hname = fs+'_invariantMass_'+era
      hists[hname] = TH1F(hname, hname, len(var_bins['inv mass']) - 1 , var_bins['inv mass'])

      #Check tf in pt bins of photon pt using 2d hist
      hname = fs+'_leadObj_pt_invariantMass_'+era
      hists[hname] = TH2F(hname, hname, len(var_bins['inv mass']) - 1 , var_bins['inv mass'], 
                                        len(var_bins['pt'])-1, var_bins['pt'])

      hname = fs+'_sublObj_pt_invariantMass_'+era
      hists[hname] = TH2F(hname, hname, len(var_bins['inv mass']) - 1 , var_bins['inv mass'], 
                                        len(var_bins['pt'])-1, var_bins['pt'])
      
      # Check tf in bins of eta
      hname = fs+'_leadObj_eta_invariantMass_'+era
      hists[hname] = TH2F(hname, hname, len(var_bins['inv mass']) - 1 , var_bins['inv mass'], 
                                        len(var_bins['eta'])-1, var_bins['eta'])
      
      hname = fs+'_sublObj_eta_invariantMass_'+era
      hists[hname] = TH2F(hname, hname, len(var_bins['inv mass']) - 1 , var_bins['inv mass'],
                                        len(var_bins['eta'])-1, var_bins['eta'])

  if isMC:
    for era in eras:
      # Count total MC events for scaling later
      hname = 'hEvents_{}'.format(era)
      hists[hname] = TH1F(hname, hname, 1, 0, 1)

      #Gen level fake rate
      for region in det_regions:
        for order in ['leadObj', 'sublObj']:
          hname= 'nEtoE_{}_{}_{}'.format(region, order, era)
          hists[hname] = TH1F(hname, hname, 1, 0, 1)
          hname= 'nEtoG_{}_{}_{}'.format(region, order, era)
          hists[hname] = TH1F(hname, hname, 1, 0, 1)

  for histName, hist in hists.items():
    hist.Sumw2()

  if subsets and isMC:
    hists = get_hists_in_subsets(dType, hists)

  return hists

def skim_egtf(filenames, dType):
  l.log(f"Running over {len(filenames)} files")

  isMC = dType != 'data'
  all_hists = make_hists(dType)

  #Store trigger indicies
  trigger_indicies = {}

  ###Start Loop over files in batch
  for filename in filenames:
    file_era = get_era(filename)

    if isMC:
      file_subset = get_subset(filename)
      hists = all_hists[file_subset]
    else:
      hists = all_hists

    l.log('Getting trigger index',1)
    trigger_name = trigger_names[file_era]
    if file_era in trigger_indicies.keys():
      trigger_index = trigger_indicies[file_era]
    else:
      trigger_index = get_trigger_index( filename, trigger_name, batch=0 if not runbatch else batch_num)
      trigger_indicies[file_era] = trigger_index

    f = TFile.Open(filename)
    t = f.Get('TreeMaker2/PreSelection')

    if isMC:
      tCounter = f.Get("tcounter")
      n_entries = tCounter.GetEntries()
      hists['hEvents_'+file_era].SetBinContent(1, n_entries)
      l.log("Total MC events: {}".format(n_entries),1)

    ##
    l.log("Starting Event Loop",1)
    n_pass = 0
    for event in t:

      eg = make_cuts(t, trigger_index, photonWP, skimmer)

      if not eg:
        continue

      n_pass+=1
      
      if isMC:
        w = t.CrossSection*t.puWeight
      else:
        w = 1

      invariantMass = (eg[0]['4vec'] + eg[1]['4vec']).M()

      fstate = ('e' if eg[0]['xseed'] else 'g') + ('B' if eg[0]['barrel'] else 'E') + \
               ('e' if eg[1]['xseed'] else 'g') + ('B' if eg[1]['barrel'] else 'E')

      hists[fstate+'_invariantMass_'+file_era].Fill(invariantMass, w)

      # Check pt dependence of photon fakes
      hists[fstate+'_leadObj_pt_invariantMass_'+file_era].Fill(invariantMass, eg[0]['pt'], w)
      hists[fstate+'_sublObj_pt_invariantMass_'+file_era].Fill(invariantMass, eg[1]['pt'], w)

      # Check eta dependence of photon fakes
      hists[fstate+'_leadObj_eta_invariantMass_'+file_era].Fill(invariantMass, eg[0]['eta'], w)
      hists[fstate+'_sublObj_eta_invariantMass_'+file_era].Fill(invariantMass, eg[1]['eta'], w)

      if isMC:
        # Gen transfer factor
        for reco, order in zip(genMatching(eg, t.GenParticles, t.GenParticles_PdgId), ['leadObj','sublObj']):
          region = 'barrel' if reco['barrel'] else 'endcap'
          if reco['genPar_name'] == 'genEle':
            if reco['xseed']:
              hists['nEtoE_{}_{}_{}'.format(region, order, file_era)].Fill(1,w)
            else:
              hists['nEtoG_{}_{}_{}'.format(region, order, file_era)].Fill(1,w)

    l.log("Closing File " + filename)
    f.Close()
    l.log(f"Number of events passing cuts: {n_pass}",1)
  return all_hists


if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(description='Skim ntuples for e-gamma transfer factor')
  parser.add_argument('dType', type=str, help='data or mc')
  parser.add_argument('--nbatches','-nb', type=int, default=1, help='number of batches')
  parser.add_argument('--batch_num','-bn', type=int, default=0, help='batch number')
  parser.add_argument('--test','-t', action='store_true', help='run in test mode')
  parser.add_argument('--verbosity','-v', type=int, default=0, help='verbosity max level')
  args = parser.parse_args()

  dType = args.dType
  batch_num = args.batch_num
  nbatches = args.nbatches
  runbatch = nbatches > 1
  test = args.test
  verbosity = args.verbosity

  isMC = dType != 'data'
  filenames = get_ntuple_filenames(dType, ntuple_version, runbatch, batch_num, nbatches)

  l.set_max_verbosity(verbosity)
  if test:
    l.set_max_verbosity(5)
    l.log("Running in test mode",1)
    # Run only one file per era
    old_filenames = filenames
    filenames = []
    for era in eras:
      era_tag = mc_eras[era] if isMC else data_eras[era]
      for filename in old_filenames:
        if era_tag in filename:
          filenames.append(filename)
          break

  hists = skim_egtf(filenames, dType)
  write_hists(hists, dType, skimmer, version, runbatch=runbatch, batch=batch_num)

