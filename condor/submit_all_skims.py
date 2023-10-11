import os

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--test', '-t', action='store_true', help='test mode')
parser.add_argument('--batch', '-b', action='store_true', help='batch mode')
args = parser.parse_args()

test = args.test

def command(hist_type, dType, nbatches):
    c = f'python submit_skims.py -s {hist_type} -d {dType} -nb {nbatches}'
    if test:
        c += ' -t'
    if args.batch:
        c += ' -b'
    return c

splitting = { 'egtf': {'data': 100,
                       'DYJetsToLL': 20},
              'closure': {'data': 100,
                          'WGJets': 20,
                          'TTGJets': 20}}

# EGTF
for hist_type, dTypes in splitting.items():
    for dType, nbatches in dTypes.items():
        if test:
            nbatches = 10000
        os.system(command(hist_type, dType, nbatches))