import htcondor
import os
import sys

import argparse

# Argparse instead
parser = argparse.ArgumentParser()
parser.add_argument('--hist_type', '-s', type=str, help='hist_type to run')
parser.add_argument('--dType', '-d', type=str, help='data or mc')
parser.add_argument('--nbatches', '-nb', type=int, default=1, help='number of batches')
parser.add_argument('--test', '-t', action='store_true', help='test mode')
parser.add_argument('--batch', '-b', action='store_true', help='batch mode')
args = parser.parse_args()

hist_type  = args.hist_type
dType    = args.dType
nbatches = args.nbatches
test     = args.test
runbatch = args.batch

#Remove condor outputs
os.system(f'rm err/{hist_type}_{dType}*.err')
os.system(f'rm log/{hist_type}_{dType}*.log')
os.system(f'rm out/{hist_type}_{dType}*.out')

def submit_batch_jobs(hist_type, dType, nbatches, batch):
  skim_events = htcondor.Submit({
      "executable": "run_skims.sh",
      "arguments": f"skim_{hist_type}.py {dType} {nbatches} {batch}",
      "output": f"out/{hist_type}_{dType}_{batch}.out",
      "error" : f"err/{hist_type}_{dType}_{batch}.err",
      "log"   : f"log/{hist_type}_{dType}_{batch}.log",              
      "request_cpus": "1",
      "request_memory": "128MB",
      "request_disk": "128MB",
  })
  print(f"Submitting batch {batch} out of {nbatches} for {hist_type}")
  schedd = htcondor.Schedd()
  submit_result = schedd.submit(skim_events)

def submit_local_jobs(hist_type, dType, nbatches, batch):
  os.system(f'./run_skims.sh skim_{hist_type}.py {dType} {nbatches} {batch}')

for batch in range(nbatches):
  if runbatch:
    submit_batch_jobs(hist_type, dType, nbatches, batch)
  else:
    submit_local_jobs(hist_type, dType, nbatches, batch)

  if test:
    break