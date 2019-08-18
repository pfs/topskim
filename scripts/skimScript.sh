#!/bin/bash

## make local directory
echo 'i am in this directory'
echo $PWD

echo ls of this directory:
ls ./


echo moving to cmssw directory ${1}

cd $1
## do cmsenv and go back
eval `scramv1 runtime -sh`

echo 'moving back to directory on node'
cd -

echo 'will run the skim config with cmsRun'
cmsRun $2


echo 'copying the output file to where it belongs'
#cp $3 $4
#use xrdcp as it's more reliable than local mounts
outD=${4/\/eos\/cms/}
xrdcp -f $3 root://eoscms//${outD}/${3}
echo 'done.'

