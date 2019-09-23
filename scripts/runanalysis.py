import os

def launch(indir,out,allowTags=[],csvWP=0):
    a= os.listdir(indir)
    for itag,tag in enumerate(a):
        if len(allowTags)>0 and not tag in allowTags: 
            continue
        if "Skim" in tag:
            extraOpts=""
        else:
            extraOpts=" --mc "
            if "Drum" in tag:
                if "amcatnlo" in tag:
                    extraOpts += " --amcatnlo "
                    if '_Muons' in tag or '_Electrons' in tag:
                        extraOpts += " --skim"
            else:
                extraOpts=" --pp "

        extraOpts += ' --csvWP %d'%csvWP
        os.system('mkdir -p skim2ljobs_csv{csv} && cd skim2ljobs_csv{csv} && python ../scripts/launchAnalysis.py {indir}/{tag} {out} {tag} {extraOpts}'.format(csv=csvWP,indir=indir, tag=tag, out=out, extraOpts=extraOpts))

for tag,csvWP in [('loose',0)]: #,('tight',1)]:
    out='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim17Sep_%s/'%tag
    indir='/eos/cms/store/cmst3/group/hintt/HIN-19-001-09Aug'
    launch(indir,out,csvWP=csvWP)

    indir='/eos/cms/store/cmst3/group/hintt/localMC'
    allowTags=['gg2ee_superchic', 'gg2mm_superchic']
    launch(indir,out,allowTags,csvWP=csvWP)
