from ROOT import TCanvas, TRatioPlot, TLegend, TH1F, TH2F, gPad, gROOT

from analysis.hists import *
from analysis.plotting import *

#Do some stuff to prove closure and things...
def bdt_closure(dType, hists):
  for era in mc_eras:
    title = '{} gg BDTclosure {}'.format(dType, era)
    plot_bdt_closure(title,
                     hists['gg_BDT_{}_{}'.format('genEG', era)],
                     hists['gg_BDT_{}_{}'.format('genGJet', era)])
    

def egtf_nonclosure_bdt(dType, tf_dict, hists):
  ''' 
  Compares gg prediction (using e->g fakerate) to actual gg in MC.
  I use this to calculate the jet nonclosure correction.
  '''
  hist_dict = { era: { bdt_bin: {} for bdt_bin in bin_names['bdt'] } for era in mc_eras}
  for era in mc_eras:
    for bdt_bin in bin_names['bdt']:
      title = '{} {} hardMET eg nonclosure {}'.format(dType, bdt_bin, era)

      hname = 'gg_{}_hardMET_{}'.format(bdt_bin, era)
      h_gg = TH1F(hname, hname, len(bins['met'])-1, bins['met'])
      hname = 'gg_prediction_{}_hardMET_{}'.format(bdt_bin, era)
      h_eg = TH1F(hname, hname, len(bins['met'])-1, bins['met'])

      #Sum over all gg object configurations
      for gg_fs, eg_fs in prediction_map.items():
        h_gg.Add( hists['{}_hardMET_{}_{}'.format(gg_fs, bdt_bin, era)])
        #Sum over all eg object configurations that contribute to this gg object configuration
        for fs in eg_fs:
          tf_region = final_state_tf_region_map[fs]    
        
          h_tmp = hists['{}_hardMET_{}_{}'.format(fs, bdt_bin, era)]
          h_tmp.Scale(tf_dict[era][tf_region])
          h_eg.Add(h_tmp)

      hist_dict[era][bdt_bin] = {'gg': h_gg, 'predicted_gg': h_eg}
      plot_eg_closure(title, h_gg, h_eg)

  return hist_dict

def egtf_closure(dType, tf_dict, hists, plot=False, save=False):
  ''' 
  Compares gg prediction (using e->g fakerate) to actual gg in MC.
  I use this to calculate the jet nonclosure correction.
  '''
  a_met_bins = analysis_bins['hardMET']
  a_bdt_bins = analysis_bins['BDT']

  h_name = lambda fs, era: '{}_hardMET_BDT_{}'.format(fs, era)

  hist_dict = { era: {} for era in mc_eras}
  for era in mc_eras:
    hname = 'gg_{}'.format(era)
    h_gg = TH2F(hname, hname, len(a_met_bins)-1, a_met_bins, len(a_bdt_bins)-1, a_bdt_bins)
    hname = 'gg_prediction_{}'.format(era)
    h_gg_pred = TH2F(hname, hname, len(a_met_bins)-1, a_met_bins, len(a_bdt_bins)-1, a_bdt_bins)

    #Sum over all gg object configurations
    for gg_fs, eg_fs in prediction_map.items():
      h_gg.Add(hists[h_name(gg_fs, era)])
      #Sum over all eg object configurations that contribute to this gg object configuration
      for fs in eg_fs:
        tf_region = final_state_tf_region_map[fs]    
        h_gg_pred.Add(hists[h_name(fs, era)], tf_dict[era][tf_region])

    hist_dict[era] = [ {'hist': h_gg, 'label': 'gg'},
                       {'hist': h_gg_pred, 'label': 'predicted gg'}]

    if plot:
      plot_eg_closure(f'{dType} gg closure'.format(era), h_gg, h_gg_pred, save=save)

  return hist_dict

def egtf_closure_genMatch(dType, tf_dict, hists):
  for era in mc_eras:
    for gg_fs, eg_fs in prediction_map.items(): #For all BDT bins
      title = '{} {} hardMET genEG closure {}'.format(dType, gg_fs, era)

      h_gg = hists['{}_hardMET_{}_{}'.format(gg_fs, 'genEG', era)]

      hname = '{}_hardMET_{}_ggPrediction_{}'.format(gg_fs, 'genEG', era)
      h_eg = TH1F(hname, hname, len(bins['met'])-1, bins['met'])
      
      for fs in eg_fs:
        tf_region = ('barrel_' if fs[fs.index('e')+1] == 'B' else 'endcap_') +\
                    ('leadFake' if fs[0] == 'e' else 'trailFake')       
        
        h_tmp = hists['{}_hardMET_{}_{}'.format(fs, 'genEG', era)]
        h_tmp.Scale(tf_dict[era][tf_region])
        h_eg.Add(h_tmp)
      r = plot_eg_closure(title, h_gg, h_eg)


    for bdt_bin in bin_names['bdt']:
      title = '{} genEG closure {} {}'.format(dType, bdt_bin, era)
      hname = 'gg_hardMET_{}_{}_{}'.format('genEG', bdt_bin, era)
      h_gg = TH1F(hname, hname, len(bins['met'])-1, bins['met'])
      hname = 'eg_hardMET_{}_prediction_{}_{}'.format('genEG', bdt_bin, era)
      h_eg = TH1F(hname, hname, len(bins['met'])-1, bins['met'])

      for gg_fs, eg_fs in prediction_map.items():
        h_gg += hists['{}_hardMET_{}_{}_{}'.format(gg_fs, 'genEG', bdt_bin, era)]
        for fs in eg_fs:
          tf_region = ('barrel_' if fs[fs.index('e')+1] == 'B' else 'endcap_') +\
                      ('leadFake' if fs[0] == 'e' else 'trailFake')       
          
          h_tmp = hists['{}_hardMET_{}_{}_{}'.format(fs, 'genEG', bdt_bin, era)]
          h_tmp.Scale(tf_dict[era][tf_region])
          h_eg += h_tmp
      plot_eg_closure(title, h_gg, h_eg)


def fake_background(dType, hists):

  h_all = hists['gg_bkg_all']
  h_fakehardMET = hists['gg_bkg_fakehardMET']
  h_realhardMET = hists['gg_bkg_realhardMET']
  h_realhardMET_lostLepton = hists['gg_bkg_realhardMET+lostLepton']

  all_bkg = h_all.Integral()

  h_all.Scale(100/all_bkg)
  h_fakehardMET.Scale(100/all_bkg)
  h_realhardMET.Scale(100/all_bkg)
  h_realhardMET_lostLepton.Scale(100/all_bkg)

  plot_background_matrix('norm_all_bkg', h_all)
  plot_background_matrix('norm_fakehardMET_bkg', h_fakehardMET)
  plot_background_matrix('norm_realhardMET_bkg', h_realhardMET)
  plot_background_matrix('norm_realhardMET_lostLepton_bkg', h_realhardMET_lostLepton)



if __name__ == '__main__': 
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('dType', type=str, help='data or mc')

  args = parser.parse_args()
  dTypes = args.dType.split(',')

  plot_dir = 'plots/closure/'

  isMC = 'data' not in dTypes
  if isMC:
    eras = mc_eras
    dType_egtf = 'DYJetsToLL'
  else:
    eras = data_eras
    dType_egtf = 'data'

  import skim_closure
  import skim_egtf
  import egtf

  if skim_closure.photonWP == skim_egtf.photonWP:
    photonWP = skim_closure.photonWP
  else:
    'PhotonWP are not the same!'
    quit()

  #Load hists
  print(eras)
  hists = skim_closure.make_hists(isMC, eras)
  for dType in dTypes:
    fn_closure = get_skim_filename(dType, 'closure')
    hists_dType = load_hists(fn_closure)
    for hName, h in hists_dType.items():
      hists[hName].Add(h)


  #Load eg transfer factor
  fn_egtf = get_skim_filename(dType_egtf, 'egtf')
  tf_dict = egtf.get_transfer_factor(fn_egtf, mc_eras[:2])
  tf_dict['Autumn18'] = tf_dict['Fall17']
  print(tf_dict)
 
  dType = dTypes[0] + '_' + dTypes[1]
  compare_bdt_fakeBKG(dType, hists)
  egtf_closure(dType, tf_dict, hists)
  egtf_nonclosure(dType, tf_dict, hists)
  fake_background(dType, hists)

