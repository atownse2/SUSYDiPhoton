#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Opyearting System on that node

# source /project01/ndcms/atownse2/SUSYDiPhoton/diphoton-env/bin/activate
source /scratch365/atownse2/SUSYDiPhoton/diphoton-env/bin/activate

cd /afs/crc.nd.edu/user/a/atownse2/Public/SUSYDiPhoton/analysis/
python skim.py $1 $2 --year $3 --nbatches $4 --batch $5 -v $6 || exit $?;
echo "job finished with status $?"