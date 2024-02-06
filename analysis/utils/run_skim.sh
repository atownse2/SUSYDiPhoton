#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Opyearting System on that node

source /scratch365/atownse2/envs/diphoton-env/bin/activate

dType=$1
analysis_region=$2
year=$3
nbatches=$4
batch=$5
verbosity=$6

PYTHON_SCRIPT_PATH=/afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/analysis/skim.py
python $PYTHON_SCRIPT_PATH $1 $2 --year $3 --nbatches $4 --batch $5 -v $6