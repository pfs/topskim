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
        os.system('mkdir -p skim1ljobs_csv{csv} && cd skim1ljobs_csv{csv} && python ../test/ttbar1l/launchAnalysis.py {indir}/{tag} {out} {tag} {extraOpts}'.format(csv=csvWP,indir=indir, tag=tag, out=out, extraOpts=extraOpts))

tag='loose'
csvWP=0
out='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim9Nov_%s/'%tag
indir='/eos/cms/store/cmst3/group/hintt/HIN-19-001-09Aug'
launch(indir,out,csvWP=csvWP)

