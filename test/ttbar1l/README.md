## Running the analysis

The selection is implemented in `bin/make1Ltree.cc`.
It uses directly as inputs the HiForest contents.
The executable can be compiled with `scram b`.
To run on a single file for testing one  can give the command
```
make1Ltree --in /eos/cms/store/cmst3/group/hintt/HIN-19-001-09Aug/TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/Chunk_0_ext0.root --out test.root --mc --amc
atnlo --max 1000
make1Ltree --in /eos/cms/store/cmst3/group/hintt/HIN-19-001-09Aug/WJetsToLNu_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_Muons/Chunk_0_ext0.root --out test_wmu.r
oot --mc --amcatnlo --max 5000
make1Ltree --in /eos/cms/store/cmst3/group/hintt/PbPb2018/TT_TuneCP5_5p02TeV-powheg-pythia8/Chunk_0_ext0.root               --out test_pp.root --pp --mc --max 1000
```
or to run on a data file:
```
make1Ltree --in /eos/cms/store/cmst3/group/hintt/HIN-19-001-09Aug/SkimMuons_04Apr2019-v1/Chunk_0_ext0.root --out test.root --max 1000
```

To loop over all the available forest trees better to use condor and then merge the outputs.
Jobs finalize in approximately 30min if queues are empty.
If needed, edit the script below for the input and output directories before running.
```
python test/ttbar1l/runanalysis.py
