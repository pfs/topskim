import ROOT
import sys
import os
import numpy as np
import pickle
from array import array
from HeavyIonsAnalysis.topskim.EventReader import *

ISOTITLES=['I','I(R=0.2)','I(R=0.25)','I(R=0.3)','Mini isolation']

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



def getVariables(dilColl):

    kin=[]
    isoEstimators=[]
    rho=[]
    for ll in dilColl:
        for i in range(2):
            l=getattr(ll,'l%d'%(i+1))
            kin.append( [l.pt,l.eta] )
            isoEstimators.append( [l.chiso+l.nhiso+l.phiso,l.isofull20,l.isofull25,l.isofull30,l.miniiso*l.pt] )
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

    qiso=np.percentile(isoEstimators,range(0,110,10),axis=0)
    qrho=np.percentile(rho,range(0,110,10))

    for i in range(5):

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
        for k in range(len(isoEstimators)):
            corIso=(isoEstimators[k][i]-isoProf.Eval(rho[k]))/kin[k][0]
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
        for k in range(len(ss_isoEstimators)):
            corIso=(ss_isoEstimators[k][i]-isoProf.Eval(ss_rho[k]))/ss_kin[k][0]
            ss_coriso.Fill(corIso)
            
        #compare both
        coriso.Scale(1./coriso.Integral())
        ss_coriso.Scale(1./ss_coriso.Integral())
        ss_coriso.SetLineColor(ROOT.kRed)
        result=getROC(coriso,ss_coriso,0.4)
        result[1].SetLineColor(1)
        result[2].SetLineColor(ROOT.kRed)
        drawIsolationProfile(hmain=coriso,
                             hextra=[(ss_coriso,'histsame')],
                             extraTxt=['#bf{[%s-UE(#rho)]/p_{T}}'%ISOTITLES[i],
                                       'I_{rel}<%f'%result[-1]],
                             name='coriso%d_%d'%(i,ch))


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kLightTemperature)        
tuneIsolation('dilepton_summary.pck',11*11)
tuneIsolation('dilepton_summary.pck',13*13)
