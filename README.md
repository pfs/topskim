# Search for PbPb->ttbar with 2018 data

## Documentation
 
See details in the analysis twiki page https://twiki.cern.ch/twiki/bin/view/CMS/PbPbTTBar2018

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

## Fitting the cross section

The physics model to be used in combine is implemented in `test/TopBtagInMediumModel.py`.
The first time it is used it should be copied to `HiggsAnalysis/CombinedLimit/python/` and compiled with `scram b`.
An example datacard is given in `test/datacard_example.txt`. 
The structure is the usual one, but notice the last lines instantiating the nuisances associated to 
the uncertainty in the b-finding parameter. 
To run the fit (creation of workspace, fit and scan of the likelihood one can do:

```
sh test/runPhysicsModelFit.sh test/datacard_example.txt
```

The likelihood scan can then be plotted with

```
sh test/drawLikelihood.sh
```


## Luminosity

```
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt
brilcalc lumi -b "STABLE BEAMS" -u /ub -i $json
```
