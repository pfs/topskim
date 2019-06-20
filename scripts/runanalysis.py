import os

def launch(indir,out,allowTags=[]):
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
            else:
                extraOpts=" --pp "
        os.system('python scripts/launchAnalysis.py {indir}/{tag} {out} {tag} {extraOpts}'.format(indir=indir, tag=tag, out=out, extraOpts=extraOpts))

out='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim19June/'

indir='/eos/cms/store/cmst3/group/hintt/HIN-19-001-v2/'
launch(indir,out)

indir='/eos/cms/store/cmst3/group/hintt/PbPb2018_rereco/'
allowTags=['SkimElectrons_04Apr2019-v1', 
           'SkimMuons_04Apr2019-v1']
launch(indir,out,allowTags)
