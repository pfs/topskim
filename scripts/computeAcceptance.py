import ROOT
import sys
import numpy as np

wgtList={'qcdscale':[1,2,3,4,6,8],
         'pdf':range(972,1075)}

fIn=ROOT.TFile.Open(sys.argv[1])
h=fIn.Get('gen_fidcounter')
h.SetDirectory(0)
hiso=fIn.Get('lep_fidcounter')
hiso.SetDirectory(0)
hfakeiso=fIn.Get('fakelep_fidcounter')
hfakeiso.SetDirectory(0)
fIn.Close()

hgen=h.ProjectionY("hgen",1,1)
hstep=[]
for xbin in xrange(1,h.GetNbinsX()+1):
    hstep.append( h.ProjectionY("hstep%d"%xbin,xbin,xbin) )
    hstep[-1].Divide(hgen)
    hstep[-1].SetDirectory(0)
    hstep[-1].Scale(100.)

    nomVal=hstep[-1].GetBinContent(1)
    print xbin,nomVal,
    for w in wgtList:
        diff=[hstep[-1].GetBinContent(i+1)-nomVal for i in wgtList[w]]
        
        if w=='qcdscale':
            print '+',max(diff),min(diff),'(scale)',
        if w=='pdf':
            print '+/-',np.std([abs(x) for x in diff]),'(pdf)',
    print ''
    

hiso.Add(hfakeiso)
hfakeiso.Divide(hiso)
fakefrac=hfakeiso.ProjectionY('fakefrac',3,3)
fakefrac.Scale(100.)

nomVal=fakefrac.GetBinContent(1)
print 'fake fraction',nomVal,
for w in wgtList:
    diff=[fakefrac.GetBinContent(i+1)-nomVal for i in wgtList[w]] 
    if w=='qcdscale':
        print '+',max(diff),min(diff),'(scale)',
    if w=='pdf':
        print '+/-',np.std([abs(x) for x in diff]),'(pdf)',
print ''
