#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node

#cd /afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/CMSSW_10_2_21/src
#eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd /afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton
source condor/diphoton-env/bin/activate

python analysis/skimmers/$1 $2 -nb $3 -bn $4
echo
echo
echo