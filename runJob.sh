#!/bin/bash
#export OUTDIR=/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis/output/v1
cd /afs/cern.ch/user/l/lvermunt/CMSSW_7_5_8_patch3/src
eval `scram r -sh`
cd /afs/cern.ch/user/l/lvermunt/CERNSummerProject/topskim
root -b -q "makeMuJetsSkim.C(\"muJetsSkim_10.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/PbPb5TeV/data/HIEWQExo/crab_FilteredSingleMuHighPt_v2/160421_135925/mergePartial/HiForest_10.root\")"
#cmsMkdir $OUTDIR
