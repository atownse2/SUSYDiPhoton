from utils import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('mc_dtypes', type=str, help='Comma separated list of MC datatypes for use in calculation of the jet normalization factor.')

args = parser.parse_args()
mc_dtypes = args.mc_dtypes.split(',')

#Get filenames
import skim_egtf
import skim_closure

if skim_closure.photonWP == skim_egtf.photonWP:
    photonWP = skim_closure.photonWP
else:
    'PhotonWP are not the same!'
    quit()

fn_closure_mc = [get_skim_filename(dType, 'closure') for dType in mc_dtypes]
fn_egtf_mc = get_skim_filename('DYJetsToLL', 'egtf')

fn_closure_data = get_skim_filename('data', 'closure')
fn_egtf_data = get_skim_filename('data', 'egtf')

###Use the MC to calculate the jet normalization factor

#Load and add MC closure hists
closure_mc_hists = skim_closure.make_hists(True, eras)
for fn in fn_closure_mc:
    hists = load_hists(fn)
    for hName, h in hists.items():
        closure_mc_hists[hName].Add(h)
closure_mc_dtype = '_'.join(mc_dtypes)

#Calculate jet normalization factor
import egtf
import closure

mc_tfs = egtf.get_transfer_factor(fn_egtf_mc, eras[:2])
mc_tfs['2018'] = mc_tfs['2017']

print('MC e->gamma transfer factors:')
print(mc_tfs)

mc_nonclosure = closure.egtf_nonclosure(closure_mc_dtype, mc_tfs, closure_mc_hists)

'''
jet_factors = { era: { bdt_bin: None for bdt_bin in bin_names['bdt'] } for era in mc_eras}

for era, bdtDict in mc_nonclosure.items():
    for bdt_bin, d in bdtDict.items():
        jet_norm = d['gg'].Integral()/d['predicted_gg'].Integral()
        jet_factors[era][bdt_bin] = jet_norm
'''

jet_factors = { era: 0 for era in eras}

for era, d in mc_nonclosure.items():
    jet_norm = d['gg'].Integral()/d['predicted_gg'].Integral()
    jet_factors[era] = jet_norm

print('Jet normalization factors:')
print(jet_factors)

### Get data prediction and scale by jet normalization factor

closure_data_hists = load_hists(fn_closure_data)
data_tfs = egtf.get_transfer_factor(fn_egtf_data, eras)

#Add all the eg final states together scaled by the appropriate transfer factors
print('Calculating data prediction...')
predictions = []
for era in eras:
    for bdt_bin in bin_names['bdt']:
        hname = 'predicted_gg_hardMET_{}_{}'.format(bdt_bin, era)
        h_eg = TH1F(hname, hname, len(bins['met'])-1, bins['met'])
        for gg_fs, eg_fs in prediction_map.items():
            for fs in eg_fs:
                tf_region = tf_region_map[fs] 
                
                h_tmp = closure_data_hists['{}_hardMET_{}_{}'.format(fs, bdt_bin, era)]
                h_tmp.Scale(data_tfs[era][tf_region])
                h_eg += h_tmp
        jet_correction = jet_factors[era]#[bdt_bin]
        h_eg.Scale(jet_correction)
        predictions.append(h_eg)

print('Saving predictions...')
### Save hists
outfile = TFile.Open('hists/ewk_predictions.root', 'recreate')
for prediction in predictions:
    overflow = drawOverflow(prediction)
    overflow.Write()
outfile.Close()
