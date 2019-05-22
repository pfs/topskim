import ROOT
import sys
import numpy as np
from optimizeLeptonIsolation import canvasHeader

def getEfficiency(fileList,flavList,ptRange,cenRange=None):

    t=ROOT.TChain('tree')
    for f in fileList : 
        t.AddFile(f)

    csv_num=ROOT.TH1F('csv_num',';CSV;Jets',50,0,1)
    csv_num.Sumw2()
    csv_den=csv_num.Clone('csvden')
    for i in range(t.GetEntries()):
        t.GetEntry(i)

        ncoll=1.
        if cenRange:
            if t.cenbin<cenRange[0] : continue
            if t.cenbin>cenRange[1] : continue

        if hasattr(t,'ncoll'):
            ncoll=t.ncoll
            if ncoll<=0 : ncoll=1.

        for j in range(t.nbjet):

            genpt=t.bjet_genpt[j]
            if genpt<ptRange[0] : continue
            if genpt>ptRange[1] : continue

            bflav=t.bjet_flavorB[j]
            if not abs(bflav) in flavList : continue

            csvVal=min(1.,max(0.,t.bjet_csvv2[j]))
            xbin=csv_num.GetXaxis().FindBin(csvVal)
            for ix in range(csv_num.GetNbinsX()+1):
                xcen=csv_num.GetXaxis().GetBinCenter(ix+1)
                if ix<xbin:
                    csv_num.Fill(xcen,ncoll)
                csv_den.Fill(xcen,ncoll)

    gr_eff=ROOT.TGraphAsymmErrors()
    gr_eff.Divide(csv_num,csv_den)
    csv_num.Delete()
    csv_den.Delete()
    return gr_eff


def showEfficiencyCurves(grColl,name):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetGridy()
    c.SetLogy()
    mg=ROOT.TMultiGraph()
    for g in grColl: mg.Add(g,'p')
    mg.Draw('ap')
    mg.GetXaxis().SetTitle('CSV')
    mg.GetYaxis().SetTitle('Efficiency')
    mg.GetXaxis().SetRangeUser(0,1)
    mg.GetYaxis().SetRangeUser(1e-4,1)
    leg=c.BuildLegend(0.15,0.15,0.4,0.3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    l=ROOT.TLine()
    l.SetLineStyle(0)
    l.SetLineColor(ROOT.kGray+1)
    for y in [0.001,0.01,0.1]:
        l.DrawLine(0,y,1,y)
    canvasHeader()
    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(name,ext))


def main():
    ROOT.gROOT.SetBatch(True)
    baseDir=sys.argv[1]
    mixSig='TTJets_TuneCP5_HydjetDrumMB-amcatnloFXFX.root'   
    ppSample='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim27Apr/TT_TuneCP5_5p02TeV-powheg-pythia8.root'
    ptRange=[30,120]

    #get the efficiency curves
    csv={}
    for tag,sample,flav,cenRange,ms,ci in [ ('b',                  baseDir+mixSig, [5],        None,     20, 1),
                                            ('b (0-30)',           baseDir+mixSig, [5],        [0,30],   20, ROOT.kRed ),
                                            ('b (30-100)',         baseDir+mixSig, [5],        [30,100], 20, ROOT.kGray ),
                                            ('b (pp)',             ppSample,       [5],        None,     20, ROOT.kGreen),
                                            ('udsg',               baseDir+mixSig, [1,2,3,21], None,     24, 1),
                                            ('udsg (0-30)',        baseDir+mixSig, [1,2,3,21], [0,30],   24, ROOT.kRed),
                                            ('udsg (30-100)',      baseDir+mixSig, [1,2,3,21], [30,100], 24, ROOT.kGray),
                                            ('udsg (pp)',          ppSample,       [1,2,3,21], None,     24, ROOT.kGreen),
                                            ('unmatched',          baseDir+mixSig, [0],        None,     21, 1),
                                            ('unmatched (0-30)',   baseDir+mixSig, [0],        [0,30],   21, ROOT.kRed),
                                            ('unmatched (30-100)', baseDir+mixSig, [0],        [30,100], 21, ROOT.kGray),
                                            ('unmatched (pp)',     ppSample,  [0],             None,     21, ROOT.kGreen)
                                            ]:
        csv[tag]=getEfficiency([sample],flav,ptRange,cenRange)
        csv[tag].SetTitle(tag)
        csv[tag].SetMarkerStyle(ms)
        csv[tag].SetLineColor(ci)
        csv[tag].SetMarkerColor(ci)

    #tune the working point
    wpEff=0.015
    bestCut=0.8
    bestEff=csv['udsg'].Eval(bestCut)
    for x in np.arange(bestCut,1,0.01):
        eff=max(csv['udsg'].Eval(x),csv['unmatched'].Eval(x))
        if abs(eff-wpEff)>abs(bestEff-wpEff): continue
        bestCut=x
        bestEff=eff
    print '<'*50
    print 'csv>',bestCut
    for tag in csv:
        print tag,'eff=',csv[tag].Eval(bestCut)
    print '<'*50

    #compare the curves
    for name,grNames in [ ('csveff',       ['b',  'b (pp)', 'udsg',  'udsg (pp)', 'unmatched', 'unmatched (pp)']),
                          ('beff',         ['b',    'b (0-30)',    'b (30-100)']),
                          ('udsgeff',      ['udsg', 'udsg (0-30)', 'udsg (30-100)']),
                          ('unmatchedeff', ['unmatched', 'unmatched (0-30)', 'unmatched (30-100)']),
                          ]:
        showEfficiencyCurves(grColl=[csv[x].Clone() for x in grNames],name=name)


if __name__ == "__main__":
    main()
