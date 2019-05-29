import ROOT
import sys
from optimizeLeptonIsolation import canvasHeader

fIn=ROOT.TFile.Open(sys.argv[1])

c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetBottomMargin(0.1)

gr={}
for p in ['e','m']:
    for d in ['pt','eta']:
        gr[(p,d)]=ROOT.TGraphAsymmErrors()
        hnum=fIn.Get('%smatch_trig_%s'%(p,d))
        hden=fIn.Get('%s_trig_%s'%(p,d))
        hnum.Rebin(4)
        hden.Rebin(4)
        gr[(p,d)].BayesDivide(hnum,hden)
        gr[(p,d)].SetName('%s_%s_trigeff'%(p,d))
        gr[(p,d)].SetMarkerStyle(20)
        gr[(p,d)].Draw('ap')
        gr[(p,d)].GetYaxis().SetRangeUser(0,1)
        canvasHeader()
        for ext in ['png','pdf']:
            c.SaveAs('%s.%s'%(gr[(p,d)].GetName(),ext))

fOut=ROOT.TFile.Open('data/trigeff_mc.root','RECREATE')
for k in gr:
    gr[k].Write()
fOut.Close()
