# Search for PbPb->ttbar with 2018 data

## Documentation
 
See details in the analysis twiki page https://twiki.cern.ch/twiki/bin/view/CMS/PbPbTTBar2018

## Running HiForrest over the grid

The full installation is as follows

```
cmsrel CMSSW_10_3_3_patch1
cd CMSSW_10_3_3_patch1/src/
cmsenv
git cms-merge-topic -u CmsHI:forest_CMSSW_10_3_1
git remote add cmshi git@github.com:CmsHI/cmssw.git
git checkout -b forest_CMSSW_10_3_1 remotes/cmshi/forest_CMSSW_10_3_1
cd HeavyIonsAnalysis/JetAnalysis/python/jets
./makeJetSequences.sh
cd -
scram build -j4
cd HeavyIonsAnalysis/JetAnalysis/test
./tests.sh
cd -
git fetch cmshi --no-tags
git checkout -b forest_CMSSW_10_3_1 remotes/cmshi/forest_CMSSW_10_3_1
git checkout forest_CMSSW_10_3_1 
git pull --rebase --no-tags
git cms-addpkg RecoHI/ZDCRecHit
git cms-addpkg RecoVertex/PrimaryVertexProducer
git cms-addpkg HLTrigger/HLTanalyzers
scram b -j 8
```

The configuration files are located in `HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_*_103X.py`
To run on the grid you can use `test/submitSkimsToCrab.py`.
You may want to edit the script for the datasets you want to process, 
the file splitting to use, output, etc. The script will produce the crab configuration files
and submit them to the crab server afterwards. 

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms
python scripts/submitSkimsToCrab.py
```


## Running the analysis

The selection/plot filling is implemented in bin/runTTto2Lselection.cc.
It uses directly as inputs the HiForest contents.
The executable can be compiled with `scram b`.
To run on a single file for testing one  can give the command
```
make2Ltree --in /eos/cms/store/cmst3/group/hintt/HIN-19-001-v2/TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/Chunk_0_ext0.root --out test.root         --mc --max 1000
make2Ltree --in /eos/cms/store/cmst3/group/hintt/PbPb2018/TT_TuneCP5_5p02TeV-powheg-pythia8/Chunk_0_ext0.root               --out test_pp.root --pp --mc --max 1000
```
or to run on a data file:
```
make2Ltree --in /eos/cms/store/cmst3/group/hintt/HIN-19-001-v2/SkimMuons_04Apr2019-v1/Chunk_0_ext0.root --out test.root --max 1000
```

To loop over all the available forest trees better to use condor and then merge the outputs.
Jobs finalize in approximately 30min if queues are empty.
If needed, edit the script below for the input and output directories before running.
```
python scripts/runanalysis.py
```

Check that all the jobs ran fine and re-run locally if needed
```
python scripts/checkCondorJobs.py ./
```

Or you can resubmit them to the batch if there are too many. this will prompt you whether you actually want to resubmit.
```
python scripts/resubmitMissing.py ./
```

Merge the outputs (hadds all the chunks according to the tag defined)
```
python scripts/mergeOutputs.py ${out}
```

Prepare data-trees for the combinatorial background from event mixing
```
python scripts/prepareCombinatorialBackgroundTree.py ${out} ${out}/Combinatorial
```

Prepare CMG tools directory structure for plotting etc.
```
python scripts/makeStructure.py ${out}
```

Optimize isolation cuts (it will use the summary pck file produced by the prepareCombinatorialBackgroundTree.py script).
```
python scripts/optimizeLeptonIsolation.py
```

Optimize b-tagging (ttbar dilepton based)
```
python scripts/optimizeBtagEff.py ${out}
```

Plot the baseline preselection plots
```
python scripts/makeAnalysisPlots.py ${out}
```

Prepare the a file with the MC trigger efficiency expectations
```
python scripts/createTrigEffSummary.py ${out}/TTJets_TuneCP5_HydjetDrumMB-amcatnloFXFX.root
```


## Luminosity

```
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt
brilcalc lumi -b "STABLE BEAMS" -u /ub -i $json
```
