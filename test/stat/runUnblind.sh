#!/bin/bash

card=$1
outdir=${card%.*}

echo making directory $outdir
echo $outdir

mkdir -p ${outdir}

echo "Creating workspace"
text2workspace.py $card -o ${outdir}/workspace.root

cd ${outdir}

commonOpts="--cminDefaultMinimizerStrategy 0 -m 172.5"
fitOpts="--robustFit=1 ${commonOpts}"
expOpts="-t -1 --expectSignal 1"

#single fit
combine -M MultiDimFit workspace.root  --algo singles --cl=0.68 ${fitOpts} -t 500 --expectSignal 1 -n expsinglestoys
combine -M MultiDimFit workspace.root  ${fitOpts} ${expOpts} --algo singles --cl=0.68 --saveFitResult -n expsingles
combine -M MultiDimFit workspace.root  ${fitOpts}            --algo singles --cl=0.68 --saveFitResult -n obssingles
combine -M FitDiagnostics workspace.root  ${fitOpts} --saveShapes --saveWithUncertainties -n obs
combine -M FitDiagnostics workspace.root  ${fitOpts} ${expOpts} --saveShapes --saveWithUncertainties -n exp

#impacts
combineTool.py -M Impacts -d workspace.root --doInitialFit        ${commonOpts} ${expOpts}
combineTool.py -M Impacts -d workspace.root --doFits --parallel 4 ${commonOpts} ${expOpts}
combineTool.py -M Impacts -d workspace.root ${commonOpts} -o impacts_exp.json
plotImpacts.py -i impacts_exp.json -o impacts_exp

combineTool.py -M Impacts -d workspace.root --doInitialFit        ${commonOpts}
combineTool.py -M Impacts -d workspace.root --doFits --parallel 4 ${commonOpts}
combineTool.py -M Impacts -d workspace.root ${commonOpts} -o impacts_obs.json
plotImpacts.py -i impacts_obs.json -o impacts_obs

#single POI likelihoods
scanOpts="--algo grid --points 300 --rMin 0 --rMax 3"
combine -M MultiDimFit workspace.root ${scanOpts} ${fitOpts} ${expOpts} -n expscan
combine -M MultiDimFit workspace.root ${scanOpts} ${fitOpts}            -n obsscan

#2D POI likelihood
scanOpts="--algo grid --points 2000 --setParameterRanges r=0,2:muzg=0.5,1.5 -P r -P muzg --fastScan"
combine -M MultiDimFit workspace.root ${scanOpts} ${fitOpts} ${expOpts} -n expscan2d  
combine -M MultiDimFit workspace.root ${scanOpts} ${fitOpts}            -n obsscan2d  

#significance
combine -M Significance workspace.root ${sigOpts} -t 1000 --expectSignal 1 -n sig_toys
combine -M Significance workspace.root ${sigOpts}                          -n sig_obs
combine -M Significance workspace.root ${expOpts}                          -n sig_exp

cd -
