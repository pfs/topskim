#!/bin/bash
cd /afs/cern.ch/user/l/lvermunt/CMSSW_7_5_8_patch3/src
eval `scram r -sh`
cd /afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/150816
#root -b -q "makeMuJetsSkim.C(\"MC_tt/MCtt_tcHigh_0.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v2/merge/HiForest_0.root\",true)"
#root -b -q "makeMuJetsSkim.C(\"MC_tt/MCtt_tcHigh_1.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v2/merge/HiForest_1.root\",true)"
root -b -q "makeMuJetsSkim.C(\"MC_tt/MCtt_tcHigh_2.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v2/merge/HiForest_2.root\",true)"
root -b -q "makeMuJetsSkim.C(\"MC_tt/MCtt_tcHigh_3.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v2/merge/HiForest_3.root\",true)"

