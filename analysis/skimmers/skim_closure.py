import sys
from array import array
import numpy as np


top = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton"
sys.path.append(top)
from analysis.selections import *
from analysis.hists import *
from analysis.metadata.sample_info import *
from analysis import logger

### Configuration

version = "101023v1"
photonWP = 'loose'

###
l = logger.Logger()
#ntuple_version = 'Skims_loosePhotonV10FullId'
ntuple_version = 'Skims_mediumPhotonV10'
hist_type = 'closure'

def make_hists(dType, subsets=True):
  '''Initialize histograms'''
  hists = OrderedDict()

  met_bins = var_bins['met']
  bdt_bins = var_bins['bdt']

  # Analysis bins
  a_met_bins = analysis_bins['hardMET']
  a_bdt_bins = analysis_bins['BDT']

  # E->Gamma Prediction and Closure Hists
  for era in eras:
    for fs in final_states:
      hname = '{}_hardMET_BDT_{}'.format(fs,era)
      hists[hname] = TH2F(hname, hname, len(a_met_bins)-1, a_met_bins, len(a_bdt_bins)-1, a_bdt_bins)

  # Count total MC events for scaling later
  if isMC:
    for era in eras:
      hname = 'hEvents_{}'.format(era)
      hists[hname] = TH1F(hname, hname, 1, 0, 1)

  # MC Closure with Object kinematics
    for era in eras:
      for fs in final_states:
        # Lead photon pT
        hname = '{}_leadPhoPt_{}'.format(fs, era)
        hists[hname] = TH1F(hname, hname, len(var_bins['pt'])-1, var_bins['pt'])
        # Sublead photon pT
        hname = '{}_sublPhoPt_{}'.format(fs, era)
        hists[hname] = TH1F(hname, hname, len(var_bins['pt'])-1, var_bins['pt'])

  # MC Closure with BDT score
    for era in eras:
      for fs in final_states:
        hname = '{}_BDT_{}'.format(fs, era)
        hists[hname] = TH1F(hname, hname, len(bdt_bins)-1, bdt_bins)

  # MC e->g Closure with Gen Matching 
    for era in eras:
      for fs in final_states:
        for genEvent in genEvents:
          hname = '{}_hardMET_{}_{}'.format(fs, genEvent, era)
          hists[hname] = TH1F(hname, hname, len(met_bins)-1, met_bins)
          for bdt_bin_name in var_bin_names['bdt']:
            hname = '{}_hardMET_{}_{}_{}'.format(fs, genEvent, bdt_bin_name, era)
            hists[hname] = TH1F(hname, hname, len(met_bins)-1, met_bins)

  # MC Comparison of fake (jet) background in BDT score
    for era in eras:
      for genEvent in genEvents:
        hname = 'gg_BDT_{}_{}'.format(genEvent, era)
        hists[hname] = TH1F(hname, hname, len(bdt_bins)-1, bdt_bins)
        for fs in gg_final_states:
          hname = '{}_BDT_{}_{}'.format(fs, genEvent, era)
          hists[hname] = TH1F(hname, hname, len(bdt_bins)-1, bdt_bins)

        # Compare BDT distributions of fake background in different MET bins
        hname = 'gg_BDT_hardMET_{}_{}'.format(genEvent, era)
        hists[hname] = TH2F(hname, hname, len(bdt_bins)-1, bdt_bins, len(met_bins)-1, met_bins)

  # Fake background matrix
    for era in eras:
      for fs in final_states:
        hname = '{}_bkg_{}'.format(fs, era)
        hists[hname] = background_matrix(hname)
        for bdt_bin_name in var_bin_names['bdt']:
          hname = '{}_bkg_{}_{}'.format(fs, bdt_bin_name, era)
          hists[hname] = background_matrix(hname)

  # gg cutflow hists with background matrix
    for cutflow in gg_cutflow:
      hname = 'gg_bkg_{}'.format(cutflow)
      hists[hname] = background_matrix(hname)

  for histName, hist in hists.items():
    hist.Sumw2()

  if subsets and isMC:
    hists = get_hists_in_subsets(dType, hists)


  return hists 
###

def skim(filenames, dType):
  '''Loop over files in filenames, fill hists'''
  l.log(f"Running over {len(filenames)} files")

  #Store trigger indicies
  trigger_indicies = {}

  all_hists = make_hists(dType)
  isMC = dType != 'data'

  for filename in filenames:
    ### Setup
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
      trigger_index = get_trigger_index(filename, trigger_name, batch=0 if not runbatch else batch_num)
      trigger_indicies[file_era] = trigger_index
    
    ### Open file
    l.log("Opening File " + filename)
    f = TFile.Open(filename)
    t = f.Get('TreeMaker2/PreSelection')

    if isMC:
      tCounter = f.Get("tcounter")
      n_entries = tCounter.GetEntries()
      hists['hEvents_'+file_era].SetBinContent(1, n_entries)
      l.log("Total MC events: {}".format(n_entries),1)

    l.log("Starting Event Loop",1)
    n_pass = 0
    for event in t:

      eg = make_cuts(t, trigger_index, photonWP, hist_type)
      if eg == False:
        continue

      n_pass += 1

      w = t.CrossSection*t.puWeight if isMC else 1
      bdt_bin = get_BDT_bin(t.mva_BDT)

      fstate = ('e' if eg[0]['xseed'] else 'g') + ('B' if eg[0]['barrel'] else 'E') + \
               ('e' if eg[1]['xseed'] else 'g') + ('B' if eg[1]['barrel'] else 'E') 
     
      MET = t.HardMETPt
      BDT = t.mva_BDT

      hists['{}_hardMET_BDT_{}'.format(fstate, file_era)].Fill(MET, BDT, w)
 
      if isMC:
        # Compare object kinematics in predicted and actual MC distributions
        hists['{}_leadPhoPt_{}'.format(fstate, file_era)].Fill(eg[0]['pt'],w)
        hists['{}_sublPhoPt_{}'.format(fstate, file_era)].Fill(eg[1]['pt'],w)
        
        # Compare BDT score in predicted and actual MC distributions
        hists['{}_BDT_{}'.format(fstate, file_era)].Fill(t.mva_BDT,w)

        #Gen Matching
        gp_names = [g['genPar_name'] for g in genMatching(eg, t.GenParticles, t.GenParticles_PdgId)]

        hists['{}_bkg_{}'.format(fstate, file_era)].Fill(gp_names[0], gp_names[1], w)
        hists['{}_bkg_{}_{}'.format(fstate, bdt_bin, file_era)].Fill(gp_names[0], gp_names[1], w)

        genEvent = False
        if 'genJet' in gp_names and 'genPho' in gp_names:
          genEvent = 'genGJet'
        elif 'genEle' in gp_names and 'genPho' in gp_names:
          genEvent = 'genEG'

        if genEvent:
          hists['{}_hardMET_{}_{}'.format(fstate, genEvent, file_era)].Fill(MET,w)
          hists['{}_hardMET_{}_{}_{}'.format(fstate, genEvent, bdt_bin, file_era)].Fill(MET,w)

        if genEvent and fstate.count('g') == 2: #gg
          hists['gg_BDT_{}_{}'.format(genEvent, file_era)].Fill(t.mva_BDT,w)
          hists['{}_BDT_{}_{}'.format(fstate, genEvent, file_era)].Fill(t.mva_BDT,w)
          hists['gg_BDT_hardMET_{}_{}'.format(genEvent, file_era)].Fill(t.mva_BDT, MET, w)
        
        if fstate.count('g') == 2: #gg
          hists['gg_bkg_{}'.format('all')].Fill(gp_names[0], gp_names[1], w)

          realhardMET = any([isNeutrino(pdgID) for pdgID in t.GenParticles_PdgId])
          if not realhardMET:
            hists['gg_bkg_{}'.format('fakehardMET')].Fill(gp_names[0], gp_names[1], w)
          else:
            hists['gg_bkg_{}'.format('realhardMET')].Fill(gp_names[0], gp_names[1], w)
            if lostLepton(t.Electrons, t.Muons, t.GenParticles, t.GenParticles_PdgId, t.GenParticles_ParentId):
              hists['gg_bkg_{}'.format('realhardMET+lostLepton')].Fill(gp_names[0], gp_names[1], w)

    # After event loop
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

  hists = skim(filenames, dType)
  write_hists(hists, dType, hist_type, version, runbatch=runbatch, batch=batch_num)
