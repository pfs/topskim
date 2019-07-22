#!/bin/bash

PO="--PO sigma_pp=69e-6 --PO A=208"
model="HiggsAnalysis.CombinedLimit.TopBtagInMediumModel:topBtags"

floatOtherPOIs="--floatOtherPOIs=1"
addOpts="--setParameters sigma=2.98 -t -1" #set to empty string to run observed

baseFitOpts="--robustFit=1 --setRobustFitStrategy=1"
baseFitOpts="${baseFitOpts} -m 172.5"
baseFitOpts="${baseFitOpts} --setParameterRanges sigma=0"
baseFitOpts="${baseFitOpts} ${floatOtherPOIs}"


echo "Make sure ${model} is in its place before running"

datacard=${1}
echo "Creating workspace from ${datacard}"

text2workspace.py ${datacard} -P ${model} ${PO} -o workspace.root

#single POI likelihoods
combine workspace.root -M MultiDimFit --algo singles --cl=0.68 ${baseFitOpts} ${addOpts} -P sigma -n sigma_singles
combine workspace.root -M MultiDimFit --algo grid --points 100 ${baseFitOpts} ${addOpts} -P sigma -n sigma_grid
