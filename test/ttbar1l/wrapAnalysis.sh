#!/bin/bash

HOME=`pwd`

#set environment
CMSSW=${1}
cd ${CMSSW}/src
eval `scram r -sh`

#run the job from the home directory
cd $HOME
input=${2}
output=${3}

extraOpts=${*:4}
opts="--in ${input} --out tto1l.root ${extraOpts}"

echo "Calling make1Ltree with [${opts}]"
make1Ltree ${opts}
cp -v tto1l.root ${output}

