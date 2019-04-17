from DataFormats.FWLite import Handle, Runs
import ROOT

lheruninfo=Handle('LHERunInfoProduct')

runs=Runs('root://cmsxrootd.fnal.gov//store/himc/RunIIpp5Spring18MiniAOD/TT_TuneCP5_5p02TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_forppRef5TeV-v1/60000/04907A40-C3B7-E811-BA32-7CD30ACDE0CC.root')
for r in runs:
    r.getByLabel('externalLHEProducer',lheruninfo)
    it=lheruninfo.product().headers_begin()
    while it!=lheruninfo.product().headers_end():
        lines=it.lines()
        allowPrint=False
        wgtCtr=0
        for i in xrange(0,lines.size()):
            linestr=lines.at(i)
            if '<weightgroup' in linestr : allowPrint=True
            if '</weightgroup' in linestr : allowPrint=False
            if not allowPrint : continue
            if 'weightgroup' in linestr :
                print '*'*50
                print linestr
                print '*'*50
            else:
                if not 'weight' in linestr : continue
                print wgtCtr,linestr
                wgtCtr+=1
        it.next()
