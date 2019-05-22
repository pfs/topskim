import os

def launch(indir,out,allowTags=[]):
    a= os.listdir(indir)
    for itag,tag in enumerate(a):
        if len(allowTags)>0 and not tag in allowTags: 
            continue
        extraOpts="true true"
        if "Skim" in tag:
            extraOpts=""
        if "Drum" in tag:
            extraOpts="true "
        os.system('python scripts/launchAnalysis.py {indir}/{tag} {out} {tag} {extraOpts}'.format(indir=indir, tag=tag, out=out, extraOpts=extraOpts))

indir='/eos/cms/store/cmst3/group/hintt/PbPb2018_rereco/'
out='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim14May/'
launch(indir,out)

indir='/eos/cms/store/cmst3/group/hintt/PbPb2018/'
allowTags=['TT_TuneCP5_5p02TeV-powheg-pythia8']
launch(indir,out,allowTags)
