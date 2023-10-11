import os
import sys

top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton"
sys.path.append(top_dir)

import json

from ROOT import TH1F, TCanvas, gROOT
import numpy as np

from collections import OrderedDict

from analysis.hists import *
from analysis.z_fit import z_fit
from analysis.metadata import sample_info

# Set TDR styles
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle();")

# Add CMS text
gROOT.LoadMacro("CMS_lumi.C")
gROOT.ProcessLine("CMS_lumi(c1);")

plot_dir = top_dir + '/plots/egtf'

def tf_BB_EE(nZ_eg, nZ_ee):
  n_ee, e_ee = nZ_ee
  n_eg, e_eg = nZ_eg

  tf = n_eg/(2*n_ee)
  e_tf = np.sqrt((e_eg/(2*n_ee))**2 + (n_eg/(2*(n_ee**2))*e_ee)**2)
  return (tf, e_tf)

def tf_BE(nZ_eg, nZ_ee):
  n_ee, e_ee = nZ_ee
  n_eg, e_eg = nZ_eg

  tf = n_eg/n_ee
  e_tf = np.sqrt((e_eg/n_ee)**2 + (n_eg/(n_ee**2)*e_ee)**2)
  return (tf, e_tf)

def tf_combined(nZ_eg, nZ_ee):
  return tf_BE(nZ_eg, nZ_ee)

def get_transfer_factors(dType, version, plot_fits=False, test=False, verbosity=0):
  '''Returns dictionary of transfer factors for all eras.

    Calculates transfer factor as a ratio of number of Z events,
    tf = nZ_egamma / nZ_ee.

    Transfer factor is derived seperately for the barrel and endcap
    regions.
    
    Additional factors of two are necessary to account for the
    fact that either electron can be faked when both objects are
    in the barrel (or endcap).
  '''
  
  # If transfer factors are saved in json, don't recalculate
  json_name = f'{top_dir}/analysis/metadata/cache/{dType}_egtf_{version}.json'
  if test: json_name = json_name.replace('.json', '_test.json')
  if (os.path.exists(json_name) and not test):
    print('Loading saved eg transfer factor')
    tf_dict = json.load(open(json_name), object_pairs_hook=OrderedDict)
    return tf_dict


  print('Calculating eg transfer factor')
  tf_dict = OrderedDict()

  fname = get_hists_filename(dType, 'egtf', version=version)
  hists = load_hists(dType, "egtf", version, verbosity=verbosity)

  eras = sample_info.eras
  if dType == "DYJetsToLL":
    eras = eras[:2]

  for era in eras:
    tf = OrderedDict()
    nZ = {}

    # BB and EE ee final states
    nZ['eBeB'] = z_fit(dType, hists['eBeB_invariantMass_{}'.format(era)], plot_fits)
    nZ['eEeE'] = z_fit(dType, hists['eEeE_invariantMass_{}'.format(era)], plot_fits)
  
    # BE ee final states
    h_eBeE_eEeB = TH1F('eBeE_eEeB_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eBeE_eEeB.Add(hists['eBeE_invariantMass_{}'.format(era)])
    h_eBeE_eEeB.Add(hists['eEeB_invariantMass_{}'.format(era)])
    nZ['eBeE_eEeB'] = z_fit(dType, h_eBeE_eEeB, plot_fits)

    # BB and EE eg final states
    h_eBgB_gBeB = TH1F('eBgB_gBeB_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eBgB_gBeB.Add(hists['eBgB_invariantMass_{}'.format(era)])
    h_eBgB_gBeB.Add(hists['gBeB_invariantMass_{}'.format(era)])
    nZ['eBgB_gBeB'] = z_fit(dType, h_eBgB_gBeB, plot_fits)

    h_eEgE_gEeE = TH1F('eEgE_gEeE_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eEgE_gEeE.Add(hists['eEgE_invariantMass_{}'.format(era)])
    h_eEgE_gEeE.Add(hists['gEeE_invariantMass_{}'.format(era)])
    nZ['eEgE_gEeE'] = z_fit(dType, h_eEgE_gEeE, plot_fits)

    # EB eg final states
    h_eEgB_gBeE = TH1F('eEgB_gBeE_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eEgB_gBeE.Add(hists['eEgB_invariantMass_{}'.format(era)])
    h_eEgB_gBeE.Add(hists['gBeE_invariantMass_{}'.format(era)])
    nZ['eEgB_gBeE'] = z_fit(dType, h_eEgB_gBeE, plot_fits)

    h_eBgE_gEeB = TH1F('eBgE_gEeB_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eBgE_gEeB.Add(hists['eBgE_invariantMass_{}'.format(era)])
    h_eBgE_gEeB.Add(hists['gEeB_invariantMass_{}'.format(era)])
    nZ['eBgE_gEeB'] = z_fit(dType, h_eBgE_gEeB, plot_fits)


    # More complicated transfer factor calculation for barrel/endcap
    h_2eBeB_eEeB_eBeE = TH1F('2eBeB_eEeB_eBeE_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_2eBeB_eEeB_eBeE.Add(hists['eBeB_invariantMass_{}'.format(era)], 2)
    h_2eBeB_eEeB_eBeE.Add(hists['eEeB_invariantMass_{}'.format(era)])
    h_2eBeB_eEeB_eBeE.Add(hists['eBeE_invariantMass_{}'.format(era)])
    nZ["2eBeB_eEeB_eBeE"] = z_fit(dType, h_2eBeB_eEeB_eBeE, plot_fits)

    h_2eEeE_eEeB_eBeE = TH1F('2eEeE_eEeB_eBeE_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_2eEeE_eEeB_eBeE.Add(hists['eEeE_invariantMass_{}'.format(era)], 2)
    h_2eEeE_eEeB_eBeE.Add(hists['eEeB_invariantMass_{}'.format(era)])
    h_2eEeE_eEeB_eBeE.Add(hists['eBeE_invariantMass_{}'.format(era)])
    nZ["2eEeE_eEeB_eBeE"] = z_fit(dType, h_2eEeE_eEeB_eBeE, plot_fits)

    h_eBgB_gBeB_gBeE_eEgB = TH1F('eBgB_gBeB_gBeE_eEgB_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eBgB_gBeB_gBeE_eEgB.Add(hists['eBgB_invariantMass_{}'.format(era)])
    h_eBgB_gBeB_gBeE_eEgB.Add(hists['gBeB_invariantMass_{}'.format(era)])
    h_eBgB_gBeB_gBeE_eEgB.Add(hists['gBeE_invariantMass_{}'.format(era)])
    h_eBgB_gBeB_gBeE_eEgB.Add(hists['eEgB_invariantMass_{}'.format(era)])
    nZ["eBgB_gBeB_gBeE_eEgB"] = z_fit(dType, h_eBgB_gBeB_gBeE_eEgB, plot_fits)

    h_eEgE_gEeE_gEeB_eBgE = TH1F('eEgE_gEeE_gEeB_eBgE_invariantMass_'+era,'Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
    h_eEgE_gEeE_gEeB_eBgE.Add(hists['eEgE_invariantMass_{}'.format(era)])
    h_eEgE_gEeE_gEeB_eBgE.Add(hists['gEeE_invariantMass_{}'.format(era)])
    h_eEgE_gEeE_gEeB_eBgE.Add(hists['gEeB_invariantMass_{}'.format(era)])
    h_eEgE_gEeE_gEeB_eBgE.Add(hists['eBgE_invariantMass_{}'.format(era)])
    nZ["eEgE_gEeE_gEeB_eBgE"] = z_fit(dType, h_eEgE_gEeE_gEeB_eBgE, plot_fits)


    # Calculate transfer factors
    tf['barrel_BB'] = tf_BB_EE(nZ['eBgB_gBeB'], nZ['eBeB'])
    tf['barrel_BE'] = tf_BE(nZ['eEgB_gBeE'], nZ['eBeE_eEeB'])
    tf['barrel_combined'] = tf_combined(nZ["eBgB_gBeB_gBeE_eEgB"], nZ["2eBeB_eEeB_eBeE"])

    tf['endcap_EE'] = tf_BB_EE(nZ['eEgE_gEeE'], nZ['eEeE'])
    tf['endcap_BE'] = tf_BE(nZ['eBgE_gEeB'], nZ['eBeE_eEeB'])
    tf['endcap_combined'] = tf_combined(nZ["eEgE_gEeE_gEeB_eBgE"], nZ["2eEeE_eEeB_eBeE"])

    # Do some rounding to make the json file more readable
    for key, values in tf.items():
      tf[key] = tuple([round(value, 4) for value in values])

    tf_dict[era] = tf

    if test:
      break

  #Save yaml
  print('Saving eg transfer factor')
  with open(json_name, 'w') as f:
    json.dump(tf_dict, f, indent=2)

  return tf_dict

def data_closure_test(fname, eras, plot_fits=False):
  # Get transfer factors
  tf_dict = get_transfer_factors(fname, eras)

  hists = load_hists(fname)

  # The number of events in the gg channels is small so add all the eras and final states together...
  h_gg = TH1F('gg_invariantMass','Invariant Mass', len(var_bins['inv mass']) - 1 , var_bins['inv mass'])
  for era in eras:
    h_gg.Add(hists['gBgB_invariantMass_{}'.format(era)])
    h_gg.Add(hists['gEgE_invariantMass_{}'.format(era)])
    h_gg.Add(hists['gBgE_invariantMass_{}'.format(era)])
    h_gg.Add(hists['gEgB_invariantMass_{}'.format(era)])
  nZ_gg = z_fit('data', h_gg, plot_fits)

  # Now add and scale the ee channels
  nZ_ee = 0
  for era in eras:
    nZ_ee += z_fit('data', hists['eBeB_invariantMass_{}'.format(era)], plot_fits)[0]*tf_dict[era]['barrel'][0]**2
    nZ_ee += z_fit('data', hists['eEeE_invariantMass_{}'.format(era)], plot_fits)[0]*tf_dict[era]['endcap'][0]**2
    nZ_ee += z_fit('data', hists['eBeE_invariantMass_{}'.format(era)], plot_fits)[0]*tf_dict[era]['barrel'][0]*tf_dict[era]['endcap'][0]
    nZ_ee += z_fit('data', hists['eEeB_invariantMass_{}'.format(era)], plot_fits)[0]*tf_dict[era]['barrel'][0]*tf_dict[era]['endcap'][0]

  # Compare
  print('Data closure test')
  print('Actual gg: {}'.format(nZ_gg))
  print('Predicted gg: {}'.format(nZ_ee))


def tf_variation(dType, version, var,
                  x_label=None, y_label=None,
                  plot_fits=False, verbosity=0):
  '''Plots transfer factor vs. var. '''
  h_name = lambda fs, var, era: f'{fs}_{var}_invariantMass_{era}'

  y_bins = None
  bin_names = None
  for b, bins in var_bins.items():
    if b in var:
      y_bins = bins
      bin_names = var_bin_names[b]
      break
  
  if y_bins is None:
    raise ValueError(f'Could not find bins for {var}')

  # Get the hists
  hists = load_hists(dType, "egtf", version)
  eras = sample_info.eras
  tf_hists = {}
  
  for era in eras:
    # Do just BB estimates for now
    h_eBeB = hists[h_name('eBeB', var, era)]
    h_eBgB_gBeB = hists[h_name('eBgB', var, era)]
    h_eBgB_gBeB.Add(hists[h_name('gBeB', var, era)])

    # Calculate transfer factor and store in TH1F
    tfs = TH1F(f'tf_{var}_{era}', 'Transfer Factor', len(y_bins) - 1, y_bins)
    for i in range(1, h_eBeB.GetNbinsY()+1):
      nZ_eBeB = z_fit(dType, h_eBeB.ProjectionX(f"eBeB_{bin_names[i]}_invariantMass_{era}",i,i,"e"), plot_fits)
      nZ_eBgB_gBeB = z_fit(dType, h_eBgB_gBeB.ProjectionX(f"eBgB_gBeB_{bin_names[i]}_invariantMass_{era}",i,i,"e"), plot_fits)
      tf = tf_BB_EE(nZ_eBgB_gBeB, nZ_eBeB)
      tfs.SetBinContent(i, tf[0])
      tfs.SetBinError(i, tf[1])
    tf_hists[era] = tfs
  
  # Plot
  plot_name = f'{dType}_egtf_{var}_{version}'

  c = TCanvas(plot_name, plot_name, 800, 600)

  y_min = 0.0
  y_max = 0.0
  for i, era in enumerate(eras):
    if era == eras[0]:
      h = tf_hists[era]
      if x_label is not None:
        h.GetXaxis().SetTitle(x_label)
      if y_label is not None:
        h.GetYaxis().SetTitle(y_label)
      h.SetLineColor(i+1)
      h.Draw('histe')
      y_min = h.GetMinimum()
      y_max = h.GetMaximum()
    else:
      tf_hists[era].SetLineColor(i+1)
      tf_hists[era].Draw('histe same')
      y_min = min(y_min, tf_hists[era].GetMinimum())
      y_max = max(y_max, tf_hists[era].GetMaximum())
  
  h.SetMinimum(y_min*0.8)
  h.SetMaximum(y_max*1.2)
  c.SaveAs(f'{plot_dir}/{plot_name}.png')
"""
def get_gen_transfer_factor(fname, eras):
  ''' Generator level fakerate for objects in this region. '''
  hists = load_hists(fname)
  fk_dict = {era:{ region:0 for region in tf_regions} for era in eras}
  for era in eras:
    for region in tf_regions:
      nEtoE = hists['nEtoE_'+region+'_'+era].GetBinContent(2) #TH1F 1 bin 0 to 1, filling with values of 1 so the overflow bin
      nEtoG = hists['nEtoG_'+region+'_'+era].GetBinContent(2)
      if nEtoG > 0:
        fk = nEtoG/(nEtoE+nEtoG)
        tf = fk/(1-fk)
        fk_dict[era][region] = tf
  return fk_dict
        
def compare_transfer_factors( tf_dict, gen_tf_dict):
  for era in eras:
    print('For era: ' + era)
    for region in tf_regions:
      print('  For region: ' + region)        
      print( '   T&P Transfer factor vs Gen Transfer factor {}:{}'.format(np.round(tf_dict[era][region], decimals=4), 
                                                                     np.round(gen_tf_dict[era][region], decimals=4)))


def plot_transfer_factor_dependence(fname, eras, var, plot_fits=False):
  ''' Plots transfer factor vs. var. '''
  dType_egtf = fname.split('_')[0]
  var_bin_names = bin_names[var] 
  hists = load_hists(fname)
  for era, eraDict in tf_dict.items():
    tf_barrel = {bin_name:0 for bin_name in var_bin_names}
    tf_endcap = {bin_name:0 for bin_name in var_bin_names}
    for bin_name in var_bin_names:
      nZ = { fs: 0 for fs in final_states}
      for fs in nZ.keys():
        h = hists[fs+'_'+bin_name+'_invariantMass_'+era]
        nZ[fs] = z_fit(dType_egtf, h, plot_fits)
      tf_b, tf_e = calculate_transfer_factor(nZ)
      tf_barrel[bin_name] = tf_b
      tf_endcap[bin_name] = tf_e
    #Plotting
    for region in ['barrel', 'endcap']:
      pass

def make_plots():
  pass
"""


if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser(description='Calculate transfer factors for eg region.')

  parser.add_argument('--dType', '-d', type=str, help='data or mc name')
  parser.add_argument('--plot_z', action='store_true', help='plot z fits')
  parser.add_argument('--test', '-t', action='store_true', help='test mode')

  args = parser.parse_args()
  dType = args.dType
  plot_z_fits = args.plot_z

  test = args.test

  version = '091323v1'

  isMC = dType != 'data'

  tf_dict = get_transfer_factors(dType, version, plot_fits=plot_z_fits)
  #data_closure_test(fn_egtf, eras, plot_fits=plot_z_fits)
  #compare_transfer_factors(tf_dict, gentf_dict)

  # if isMC:
    # gentf_dict = get_gen_transfer_factor(fn_egtf, eras)


