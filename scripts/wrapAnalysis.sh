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
localout=`basename ${output}`

extraOpts=""
isMC=${4}
if [ ! -z "${isMC}" ]; then
   extraOpts="${extraOpts} --mc"
fi
isPP=${5}
if [ ! -z "${isPP}" ]; then
    extraOpts="${extraOpts} --pp"
fi
opts="--in ${input} --out ${localout} ${extraOpts}"

echo "Calling runTTto2Lselection with [${opts}]"
runTTto2Lselection ${opts}
cp -v ${localout} ${output}

