import os

indir='/eos/cms/store/cmst3/group/hintt/PbPb2018_rereco/'
a= os.listdir(indir)
out='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim14May/'
for itag,tag in enumerate(a):
    extraOpts="true true"
    if "Skim" in tag:
       extraOpts=""
    if "Drum" in tag:
        extraOpts="true "
    os.system('python scripts/launchAnalysis.py {indir}/{tag} {out} {tag} {extraOpts}'.format(indir=indir, tag=tag, out=out, extraOpts=extraOpts))
