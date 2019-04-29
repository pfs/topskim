# Search for PbPb->ttbar with 2018 data

## Documentation
 
See details in the analysis twiki page https://twiki.cern.ch/twiki/bin/view/CMS/PbPbTTBar2018

## Running the analysis

The selection/plot filling is implemented in bin/runTTto2Lselection.cc.
It uses directly as inputs the HiForest contents.
The executable can be compiled with `scram b`.
To run on a single file for testing one  can give the command
```
make2Ltree --in /eos/cms/store/cmst3/group/hintt/PbPb2018/TTDilepton_TuneZ2_HydjetDrumMB/Chunk_0_ext0.root --out test.root --mc --max 1000
```
or to run on a data file:
```
make2Ltree --in /eos/cms/store/cmst3/group/hintt/PbPb2018/SkimMuons_PromptRecov1/Chunk_0_ext0.root --out test.root --max 1000
```

To loop over all the available forest trees better to use condor and then merge the outputs.
Jobs finalize in approximately 30min if queues are empty.
```
dir=/eos/cms/store/cmst3/group/hintt/PbPb2018/
a=(`ls ${dir}`)
out=/eos/cms/store/cmst3/group/hintt/PbPb2018_skim27Apr
for tag in ${a[@]}; do    
    extraOpts="true true"
    if [[ $tag == *"Skim"* ]]; then
       extraOpts=""
    fi
    if [[ $tag == *"Drum"* ]]; then
       extraOpts="true"
    fi
    python scripts/launchAnalysis.py ${dir}/${tag} $out $tag $extraOpts;
done   
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

Optimize isolation cuts (it will use the summary pck file produced by the prepareCombinatorialBackgroundTree.py script).
```
python scripts/optimizeLeptonIsolation.py ${out}
```

Plot the baseline preselection plots
```
python scripts/makeAnalysisPlots.py plots/
```

## Luminosity

```
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt
brilcalc lumi -b "STABLE BEAMS" -u /ub -i $json
```
