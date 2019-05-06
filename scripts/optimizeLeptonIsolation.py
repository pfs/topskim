import ROOT
import sys
import os
import numpy as np
import pickle
from array import array
from HeavyIonsAnalysis.topskim.EventReader import *
from prepareCombinatorialBackgroundTree import prepareDileptonCollection

ISOTITLES=['I_{raw}','I','I(R=0.2)','I(R=0.25)','I(R=0.3)','Mini isolation','Charged isolation']
REQ_EFFB=0.4
MCTAG='TTDilepton_TuneZ2_HydjetDrumMB'

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

def tuneIsolation(mixFile,ch,matchedll=None):
    
    #readout Z and same-sign candidate events
    with open(mixFile,'r') as cache:
        allDileptons=pickle.load(cache)
        zll=allDileptons[(ch,True)]
        ssll=[ ll for ll in allDileptons[(ch,False)] if ll.isSS ] 

    kin,isoEstimators,rho=getVariables(zll)
    ss_kin,ss_isoEstimators,ss_rho=getVariables(ssll)
    mc_kin,mc_isoEstimators,mc_rho=None,None,None
    if matchedll:
        mc_kin,mc_isoEstimators,mc_rho=getVariables(matchedll)

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
            if i>0:
                corIso=(isoEstimators[k][i]-isoProf.Eval(rho[k]))/kin[k][0]
            else:
                corIso=(isoEstimators[k][i])/kin[k][0]
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
            if i>0:
                corIso=(ss_isoEstimators[k][i]-isoProf.Eval(ss_rho[k]))/ss_kin[k][0]
            else:
                corIso=(ss_isoEstimators[k][i])/ss_kin[k][0]
            ss_corIsoEstimators.append(corIso)
            ss_coriso.Fill(corIso)
            

        #repeat the same for MC truth if available
        mc_coriso=None
        mc_corIsoEstimators=[]
        if mc_isoEstimators:
            mc_coriso=coriso.Clone('mccoriso%d'%i)
            mc_coriso.Reset('ICE')
            mc_coriso.SetLineWidth(2)
            mc_coriso.SetLineColor(ROOT.kCyan+1)
            mc_coriso.SetLineStyle(9)
            for k in range(len(mc_isoEstimators)):
                if i>0:
                    corIso=(mc_isoEstimators[k][i]-isoProf.Eval(mc_rho[k]))/mc_kin[k][0]
                else:
                    corIso=(mc_isoEstimators[k][i])/mc_kin[k][0]
                mc_corIsoEstimators.append(corIso)
                mc_coriso.Fill(corIso)
            mc_coriso.Scale(1./mc_coriso.Integral())

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
        hextra=[(ss_coriso,'histsame')]
        if mc_coriso: hextra.append( (mc_coriso,'histsame') )
        drawIsolationProfile(hmain=coriso,
                             hextra=hextra,
                             extraTxt=[isoTitle,'I_{rel}<%f'%cut],
                             name='coriso%d_%d'%(i,ch))

        #efficiency versus lepton kinematics
        isoeff=ROOT.TH2F('effvsptvseta_iso%d'%i, ';Transverse momentum [GeV];Pseudo-rapidity;Efficiency', 5,array('d',qkin[:,0]),5,array('d',qkin[:,1]) )
        isoeff.Sumw2()
        isoeff_den=isoeff.Clone('isoeffden')    
        for k in range(len(mc_kin)):
            pt,eta,_=mc_kin[k]
            isoeff_den.Fill(pt,eta)
            if corIsoEstimators[k]>cut : continue
            isoeff.Fill(pt,eta)
        isoeff.Divide(isoeff_den)
        isoeff_den.Delete()
        drawIsolationProfile(hmain=isoeff,
                             hextra=[],
                             extraTxt=[isoTitle,'I_{rel}<%f'%cut],                                       
                             name='iso%deff_%d'%(i,ch))

        #scale factor
        if mc_isoEstimators:
            mcisoeff=isoeff.Clone('isomcnum')
            mcisoeff.Reset('ICE')
            mcisoeff_den=mcisoeff.Clone('isomcden')
            for k in range(len(mc_kin)):
                pt,eta,_=mc_kin[k]
                mcisoeff_den.Fill(pt,eta)
                if mc_corIsoEstimators[k]>cut : continue
                mcisoeff.Fill(pt,eta)
            mcisoeff.Divide(mcisoeff_den)
            mcisoeff_den.Delete()
            sfisoeff=isoeff.Clone('sfiso')
            sfisoeff.Divide(mcisoeff)
            sfisoeff.GetZaxis().SetTitle('SF = #varepsilon(data)/#varepsilon(MC)')
            #sfisoeff.GetZaxis().SetRangeUser(0.8,1.2)
            drawIsolationProfile(hmain=sfisoeff,
                             hextra=[],
                             extraTxt=[isoTitle,'I_{rel}<%f'%cut],                                       
                             name='sfiso%deff_%d'%(i,ch))

            
        #look for hotspots in same-sign data
        ssHotSpots=ROOT.TH2F('ssetavsphi',';Pseudo-rapidity;Azimuthal angle [rad];Candidates',25,-2.5,2.5,25,-3.15,3.15)
        isossHotSpots=ssHotSpots.Clone('isossetavsphi')
        if i==5:
            sshs_iso=ROOT.TH1F('sshsiso%d'%i,';%s;Leptons'%ISOTITLES[i],10,array('d',qiso[:,i]))
        else:
            sshs_iso=ROOT.TH1F('sshsiso%d'%i,';%s;Leptons'%ISOTITLES[i],50,-3,3)
        ssnhs_iso=sshs_iso.Clone('ssnhsiso%d'%i)
        ziso=sshs_iso.Clone('ziso%d'%i)
        for k in range(len(ss_kin)):
            
            pt,eta,phi=ss_kin[k]
            HEM1516failure=True if eta<-1.39 and -1.6<phi<-0.9 else False

            isoVal=ss_isoEstimators[k][i] if i==5 else ss_corIsoEstimators[k]
            if HEM1516failure:
                sshs_iso.Fill(isoVal)
            else:
                ssnhs_iso.Fill(isoVal)

            ssHotSpots.Fill(eta,phi)
            if ss_corIsoEstimators[k]>cut : continue
            isossHotSpots.Fill(eta,phi)

        for k in range(len(kin)):
            isoVal=isoEstimators[k][i] if i==5 else corIsoEstimators[k]
            ziso.Fill(isoVal)


        if i==0: drawIsolationProfile(hmain=ssHotSpots,hextra=[],extraTxt=[],name='sshotspots_%d'%ch)
        drawIsolationProfile(hmain=isossHotSpots,hextra=[],extraTxt=[isoTitle,'I_{rel}<%f'%cut],name='iso%dsshotspots_%d'%(i,ch))
        ssHotSpots.Delete()
        isossHotSpots.Delete()

        if ch==121:
            for h in [ssnhs_iso,sshs_iso,ziso]:
                h.Scale(1./h.Integral())
                h.SetLineWidth(2)
                h.GetYaxis().SetTitle('PDF')
                if i!=5: continue
                for xbin in range(ssnhs_iso.GetNbinsX()):
                    wid=ssnhs_iso.GetXaxis().GetBinWidth(xbin+1)                
                    h.SetBinContent(xbin+1,h.GetBinContent(xbin+1)/wid)
                    h.SetBinError(xbin+1,h.GetBinError(xbin+1)/wid)
            ziso.SetLineColor(ROOT.kRed)
            sshs_iso.SetLineColor(ROOT.kGray+2)
            sshs_iso.SetFillStyle(3001)
            sshs_iso.SetFillColor(ROOT.kGray)
            drawIsolationProfile(hmain=ziso,
                                 hextra=[(sshs_iso,'histsame'),(ssnhs_iso,'histsame')],
                                 extraTxt=[ISOTITLES[i] if i==5 else isoTitle],
                                 name='iso%d_hotspotcomp_%d'%(i,ch))

    #overall performance summary
    drawROCSummary(rocCurves,'rocsummary_%d'%ch)

def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)        

    matchedll={11*11:None,13*13:None}
    if len(sys.argv)>1:
       
        print 'Processing MC truth from',MCTAG
        #prepareDileptonCollection(sys.argv[1],MCTAG)

        #readout matched leptons
        with open('dilepton_summary_%s.pck'%MCTAG,'r') as cache:
            allDileptons=pickle.load(cache)
            for ch in matchedll:
                matchedll[ch]=[ll for ll in allDileptons[(ch,False)] if ll.l1.matched and ll.l2.matched]

    for ch in matchedll:
        tuneIsolation('dilepton_summary.pck',ch,matchedll[ch])

if __name__ == "__main__":
    main()
