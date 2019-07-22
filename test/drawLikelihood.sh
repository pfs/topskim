#!/bin/bash

python scripts/plotLikelihoodScanResults.py \
    "Expected"=higgsCombinesigma_grid.MultiDimFit.mH172.5.root \
    -p sigma;


#python plotLikelihoodScanResults.py \
#    "Expected"=higgsCombinesigma_eb_grid.MultiDimFit.mH172.5.root \
#    -p sigma,eb

