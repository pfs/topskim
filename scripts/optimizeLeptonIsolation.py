import ROOT
import sys
import os
import numpy as np
import pickle
from array import array
from HeavyIonsAnalysis.topskim.EventReader import *

ISOTITLES=['I','I(R=0.2)','I(R=0.25)','I(R=0.3)','Mini isolation','Charged isolation']
REQ_EFFB=0.4

def canvasHeader(extraTxt=[]):
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    txt.DrawLatex(0.95,0.97,'#scale[0.8]{1618 #mub^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')
    txt.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
    txt.SetTextSize(0.035)
    for i in range(len(extraTxt)):
        txt.DrawLatex(0.15,0.9-i*0.04,extraTxt[i])


def drawIsolationProfile(hmain,hextra=[],extraTxt=[],name='isoprof'):

    """ draws the isolation profile """
   
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    if hmain.InheritsFrom('TH2'):
        c.SetRightMargin(0.12)
        hmain.Draw('colz')
    elif hmain.InheritsFrom('TH1'):
        c.SetRightMargin(0.03)
        hmain.Draw('hist')
    for h,drawOpt in hextra:
        h.Draw(drawOpt)
    canvasHeader(extraTxt)
    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']:
        c.SaveAs('{0}.{1}'.format(name,ext))


def drawROCSummary(grColl,name):

    colors=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]
    ncols=len(colors)
    lines=[1,2,3]

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetGridy()

    mg=ROOT.TMultiGraph()
    for i in range(len(grColl)):
        grColl[i].SetLineColor( colors[i%ncols] )
        grColl[i].SetLineStyle( lines[i/ncols] )
        grColl[i].SetLineWidth(2)
        mg.Add(grColl[i],'l')
    mg.Draw('al')
    mg.GetXaxis().SetTitle('Signal efficiency')
    mg.GetYaxis().SetTitle('Background efficiency')
    mg.GetXaxis().SetRangeUser(0.5,1)
    mg.GetYaxis().SetRangeUser(0.2,1)
    leg=c.BuildLegend(0.15,0.94,0.4,0.94-0.03*len(grColl))
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)   
    l=ROOT.TLine(0.5,REQ_EFFB,1,REQ_EFFB)
    l.SetLineColor(ROOT.kGray+2)
    l.Draw()    
    canvasHeader()
    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']:
        c.SaveAs("%s.%s"%(name,ext))


def getVariables(dilColl):

    kin=[]
    isoEstimators=[]
    rho=[]
    for ll in dilColl:
        for i in range(2):
            l=getattr(ll,'l%d'%(i+1))
            kin.append( [l.pt,l.eta,l.phi] )
            isoEstimators.append( [l.chiso+l.nhiso+l.phiso,l.isofull20,l.isofull25,l.isofull30,l.miniiso*l.pt,l.chiso] )
            rho.append( l.rho )

    return kin,isoEstimators,rho


def getROC(sig, bkg, cutAtEffB):

    """build the efficiency and ROC curves"""

    #start the graphs
    roc_gr=ROOT.TGraph()
    roc_gr.SetName('roc')
    effs_gr=roc_gr.Clone('effs')
    effb_gr=roc_gr.Clone('effb')

    #loop over bins
    tot_sig=sig.Integral()
    tot_bkg=bkg.Integral()
    nbins=sig.GetNbinsX()
    bestEff,bestCut=None,None
    for xbin in range(nbins+1):
        cut=sig.GetXaxis().GetBinUpEdge(xbin)
        
        effsig=sig.Integral(0,xbin+1)/tot_sig if tot_sig>0 else 0
        effbkg=bkg.Integral(0,xbin+1)/tot_bkg if tot_bkg>0 else 0

        #if not bestEff or abs(effsig-cutAtEffS)<abs(bestEff-cutAtEffS):
        #    bestEff=effsig
        #    bestCut=cut

        if not bestEff or abs(effbkg-cutAtEffB)<abs(bestEff-cutAtEffB):
            bestEff=effbkg
            bestCut=cut

        roc_gr.SetPoint(xbin-1,effsig,effbkg)
        effs_gr.SetPoint(xbin-1,cut,effsig)
        effb_gr.SetPoint(xbin-1,cut,effbkg)

    return roc_gr,effs_gr,effb_gr, bestEff, bestCut

def tuneIsolation(mixFile,ch=11*11):
    
    #readout Z and same-sign candidate events
    with open(mixFile,'r') as cache:
        allDileptons=pickle.load(cache)
        zll=allDileptons[(ch,True)]
        ssll=[ ll for ll in allDileptons[(ch,False)] if ll.isSS ] 

    kin,isoEstimators,rho=getVariables(zll)
    ss_kin,ss_isoEstimators,ss_rho=getVariables(ssll)

    qkin=np.percentile(kin,range(0,120,20),axis=0)
    qiso=np.percentile(isoEstimators,range(0,110,10),axis=0)
    qrho=np.percentile(rho,range(0,110,10))

    rocCurves=[]
    for i in range(6):

        #parametrize the profile versus rho
        isoProf=ROOT.TF1('isoprof','[0]*pow(x+[1],2)+[2]*(x+[1])',0,1000)
        isovsrho=ROOT.TH2F('iso%dvsrho'%i, ';#rho;%s [GeV]'%ISOTITLES[i], 10,array('d',qrho),10,array('d',qiso[:,i]) )
        for k in range(len(isoEstimators)):
            isovsrho.Fill(rho[k],isoEstimators[k][i])
        isovsrho_prof=isovsrho.ProfileX()
        isovsrho_prof.Fit(isoProf,'MRQ+')
        isovsrho_prof.SetMarkerStyle(20)
        drawIsolationProfile(hmain=isovsrho,
                             hextra=[(isovsrho_prof,'e1same')],
                             extraTxt=['#bf{%s}'%ISOTITLES[i],
                                       'UE=a(#rho+b)^{2}+c(#rho+b)',
                                       'a=%f'%isoProf.GetParameter(0),
                                       'b=%f'%isoProf.GetParameter(1),
                                       'c=%f'%isoProf.GetParameter(2)],
                             name='iso%dprof_%d'%(i,ch))


        #correct the isolation estimator
        corisovsrho=ROOT.TH2F('coriso%dvsrho'%i, ';#rho;[%s-UE(#rho)]/p_{T}'%ISOTITLES[i], 10,array('d',qrho),50,-3,3)
        coriso=ROOT.TH1F('coriso%d'%i, ';[%s-UE(#rho)]/p_{T};'%ISOTITLES[i],50,-3,3)
        corIsoEstimators=[]
        for k in range(len(isoEstimators)):
            corIso=(isoEstimators[k][i]-isoProf.Eval(rho[k]))/kin[k][0]
            corIsoEstimators.append(corIso)
            corisovsrho.Fill(rho[k],corIso)
            coriso.Fill(corIso)
        corisovsrho_prof=corisovsrho.ProfileX()
        corisovsrho_prof.SetMarkerStyle(20)
        drawIsolationProfile(hmain=corisovsrho,
                             hextra=[(corisovsrho_prof,'e1same')],
                             extraTxt=['#bf{[%s-UE(#rho)]/p_{T}}'%ISOTITLES[i]],
                             name='coriso%dprof_%d'%(i,ch))

        #repeat for the same-sign events
        ss_coriso=coriso.Clone('sscoriso%d'%i)
        ss_coriso.Reset('ICE')
        ss_corIsoEstimators=[]
        for k in range(len(ss_isoEstimators)):
            corIso=(ss_isoEstimators[k][i]-isoProf.Eval(ss_rho[k]))/ss_kin[k][0]
            ss_corIsoEstimators.append(corIso)
            ss_coriso.Fill(corIso)
            
        #compare both
        coriso.Scale(1./coriso.Integral())
        coriso.SetLineWidth(2)
        coriso.GetYaxis().SetTitle('PDF')
        ss_coriso.Scale(1./ss_coriso.Integral())
        ss_coriso.SetLineColor(ROOT.kGray+2)
        ss_coriso.SetFillStyle(3001)
        ss_coriso.SetFillColor(ROOT.kGray)
        result=getROC(coriso,ss_coriso,REQ_EFFB)
        rocCurves.append(result[0])
        isoTitle='#bf{[%s-UE(#rho)]/p_{T}}'%ISOTITLES[i]
        rocCurves[-1].SetTitle(isoTitle)
        result[1].SetLineColor(1)
        result[2].SetLineColor(ROOT.kRed)
        cut=result[-1]
        drawIsolationProfile(hmain=coriso,
                             hextra=[(ss_coriso,'histsame')],
                             extraTxt=[isoTitle,'I_{rel}<%f'%cut],
                             name='coriso%d_%d'%(i,ch))

        #efficiency versus lepton kinematics
        isoeff=ROOT.TH2F('effvsptvseta_iso%d'%i, ';Transverse momentum [GeV];Pseudo-rapidity;Efficiency', 5,array('d',qkin[:,0]),5,array('d',qkin[:,1]) )
        isoeff.Sumw2()
        isoeff_den=isoeff.Clone('isoeffden')
        for k in range(len(kin)):
            isoeff_den.Fill(kin[k][0],kin[k][1])
            if corIsoEstimators[k]>cut : continue
            isoeff.Fill(kin[k][0],kin[k][1])
        isoeff.Divide(isoeff_den)
        isoeff_den.Delete()
        drawIsolationProfile(hmain=isoeff,
                             hextra=[],
                             extraTxt=[isoTitle,'I_{rel}<%f'%cut],                                       
                             name='iso%deff_%d'%(i,ch))


        #look for hotspots in same-sign data
        ssHotSpots=ROOT.TH2F('ssetavsphi',';Pseudo-rapidity;Azimuthal angle [rad];Candidates',25,-2.5,2.5,25,-3.15,3.15)
        isossHotSpots=ssHotSpots.Clone('isossetavsphi')
        for k in range(len(ss_kin)):
            ssHotSpots.Fill(ss_kin[k][1],ss_kin[k][2])
            if ss_corIsoEstimators[k]>cut : continue
            isossHotSpots.Fill(ss_kin[k][1],ss_kin[k][2])
        if i==0: drawIsolationProfile(hmain=ssHotSpots,hextra=[],extraTxt=[],name='sshotspots_%d'%ch)
        drawIsolationProfile(hmain=isossHotSpots,hextra=[],extraTxt=[isoTitle,'I_{rel}<%f'%cut],name='iso%dsshotspots_%d'%(i,ch))
        ssHotSpots.Delete()
        isossHotSpots.Delete()

    #overall performance summary
    drawROCSummary(rocCurves,'rocsummary_%d'%ch)

def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)        
    tuneIsolation('dilepton_summary.pck',11*11)
    tuneIsolation('dilepton_summary.pck',13*13)

if __name__ == "__main__":
    main()
