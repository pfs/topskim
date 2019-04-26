import ROOT
import sys
import os
import numpy as np
import pickle
from HeavyIonsAnalysis.topskim.EventReader import *


def tuneIsolation(mixFile,ch=11*11):
    
    #readout Z candidate events
    zll=None
    ss=None
    with open(mixFile,'r') as cache:
        zll=pickle.load(cache)[ch,True]
    
    isoEstimators=[]
    for ll in zll:
        for i in range(2):
            l=getattr(ll,'l%d'%(i+1))
            isoEstimators.append( [l.chiso+l.nhiso+l.phiso,l.isofull20,l.isofull25,l.isofull30,l.miniiso*l.pt,l.rho] )

    q=np.percentile(isoEstimators,range(0,110,10),axis=0)

    histos={}
    for i in range(5):
        histos['isovsrho%d'%i]=ROOT.TH2F('%s_%svs%s'%(isoKey,iso,rho), ';%s;%s'%(RHOTITLE[rho],ISOTITLE[iso]), 10,array('d',q[:,1]),10,array('d',q[:,0]))



tuneIsolation('dilepton_summary.pck')
#
#    
#    #start the total isolation estimator per lepton
#    isoVals=[0]*len(leptons)
#
#    #determine the isolation components separately
#    histos={}
#    for iso, rho in isoComponents:
#
#        isoCompVals=[]
#        for l,isZ,isSS,cenbin in leptons:
#            if not isZ: continue
#            iiso=getattr(l,iso)
#            irho=l.rho[rho] if rho else 0.
#            isoCompVals.append([iiso,irho,cenbin])
#
#        #parameterize rho vs centrality and isolation versus rho
#        if rho:
#
#            q=np.percentile(isoCompVals,range(0,110,10),axis=0)
#            histos[iso]={
#                'rhovscen':ROOT.TH2F('%s_%svscen'%(isoKey,rho), ';Centrality;%s'%RHOTITLE[rho], 5,array('d',zcenq), 10,array('d',q[:,1])),
#                'isovsrho':ROOT.TH2F('%s_%svs%s'%(isoKey,iso,rho), ';%s;%s'%(RHOTITLE[rho],ISOTITLE[iso]), 10,array('d',q[:,1]),10,array('d',q[:,0]))
#                }
#            for iiso,irho,cenbin in isoCompVals: 
#                histos[iso]['rhovscen'].Fill(cenbin,irho)
#                histos[iso]['isovsrho'].Fill(irho,iiso)
#
#            histos[iso]['rhovscen_prof']=histos[iso]['rhovscen'].ProfileX()
#            histos[iso]['rhovscen_prof'].Fit(rhoProf,'MRQ+')
#            histos[iso]['rhovscen_prof_form']='#rho=%.g/(%.g+c)'%(rhoProf.GetParameter(0),rhoProf.GetParameter(1))
#
#            histos[iso]['isovsrho_prof']=histos[iso]['isovsrho'].ProfileX()
#            histos[iso]['isovsrho_prof'].Fit(isoProf,'MRQ+')
#            histos[iso]['isovsrho_prof_form']='UE=%f(#rho+%f)^{2}+%f(#rho+%f)'%(isoProf.GetParameter(0),isoProf.GetParameter(1),isoProf.GetParameter(2),isoProf.GetParameter(1))
#            
#        #compute isolation component
#        for i in range(len(leptons)):
#            iiso=getattr(leptons[i][0],iso)
#            if rho: iiso-=isoProf.Eval(leptons[i][0].rho[rho])
#            isoVals[i]+=iiso 
#
#
#    #profile the absolute isolation versus centrality (signal and background)
#    q=np.percentile(isoVals,range(0,110,10),axis=0)
#    histos[isoKey]={
#        'isovscen'   : ROOT.TH2F('%svscen'%isoKey,'Z#rightarrowll;Centrality;%s [GeV]'%ISOTITLE[isoKey], 5,array('d',zcenq), 10,array('d',q)),
#        'bkgisovscen': ROOT.TH2F('bkg%svscen'%isoKey,'same-sign ll;Centrality;%s [GeV]'%ISOTITLE[isoKey], 5,array('d',zcenq), 10,array('d',q)),
#        }
#    for i in range(len(leptons)):
#        l,isZ,isSS,cenbin=leptons[i]
#        histos[isoKey]['isovscen' if isZ else 'bkgisovscen'].Fill(cenbin,isoVals[i])
#    histos[isoKey]['isovscen_prof']=histos[isoKey]['isovscen'].ProfileX()
#    histos[isoKey]['isovscen_prof'].Fit(cenProf,'MRQ+')
#    histos[isoKey]['isovscen_prof_form']='I=exp^{%.gc}[%.gc^{2}+%.gc+%.g]'%(cenProf.GetParameter(0),cenProf.GetParameter(1),cenProf.GetParameter(2),cenProf.GetParameter(3))
#
#    #final plots for relative isolation with and without centrality correction
#    riRan=(-5,5)
#    #if 'miniiso' in isoKey : riRan=(-1,5)
#    histos['rel'+isoKey]={
#        'isovscen_rel'       : ROOT.TH2F('%svscen_rel'%isoKey,        'Z#rightarrowll;Centrality;%s/p_{T}'%ISOTITLE[isoKey], 5,array('d',zcenq), 200,riRan[0],riRan[1]),
#        'bkgisovscen_rel'    : ROOT.TH2F('bkg%svscen_rel'%isoKey,     'same-sign ll;Centrality;%s/p_{T}'%ISOTITLE[isoKey],   5,array('d',zcenq), 200,riRan[0],riRan[1]),
#        'isovscen_relccor'   : ROOT.TH2F('%svscen_relccor'%isoKey,    'Z#rightarrowll;Centrality;%s/p_{T}'%ISOTITLE[isoKey], 5,array('d',zcenq), 200,riRan[0],riRan[1]),
#        'bkgisovscen_relccor': ROOT.TH2F('bkg%svscen_relccor'%isoKey, 'same-sign ll;Centrality;%s/p_{T}'%ISOTITLE[isoKey],   5,array('d',zcenq), 200,riRan[0],riRan[1])
#        }
#    for i in range(len(leptons)):
#        l,isZ,isSS,cenbin=leptons[i]
#        pt=l.p4.Pt()
#        iso=isoVals[i]
#        corIso=iso-cenProf.Eval(cenbin)
#        histos['rel'+isoKey]['%sisovscen_rel'%('' if isZ else 'bkg')].Fill(cenbin,iso/pt)
#        histos['rel'+isoKey]['%sisovscen_relccor'%('' if isZ else 'bkg')].Fill(cenbin,corIso/pt)
#
#    
#    effCurves={'' :getROC(histos['rel'+isoKey]['isovscen_rel'].ProjectionY('sig'),
#                          histos['rel'+isoKey]['bkgisovscen_rel'].ProjectionY('bkg'),
#                          isoKey, ISOTITLE[isoKey], ISOCOLOR[isoKey], 1),
#               'ccor':getROC(histos['rel'+isoKey]['isovscen_relccor'].ProjectionY('sig'),
#                             histos['rel'+isoKey]['bkgisovscen_relccor'].ProjectionY('bkg'),
#                             isoKey,
#                             ISOTITLE[isoKey]+'*', ISOCOLOR[isoKey],2) }
#
#    effCenCurves={}
#    for tag in ['','ccor']:
#
#        cut=effCurves[tag][-1]
#        effCenCurves[tag]=(effCurves[tag][1].Clone('sigvscenprof'+isoKey),
#                           effCurves[tag][2].Clone('bkgvscenprof'+isoKey))
#        effCenCurves[tag][0].Set(0)
#        effCenCurves[tag][1].Set(0)
#
#        nbinsy=histos['rel'+isoKey]['isovscen_rel'+tag].GetNbinsY()
#        nbinsx=histos['rel'+isoKey]['isovscen_rel'+tag].GetNbinsX()
#        for xbin in range(nbinsx):
#            cen=histos['rel'+isoKey]['isovscen_rel'+tag].GetXaxis().GetBinCenter(xbin+1)
#
#            ybin=histos['rel'+isoKey]['isovscen_rel'+tag].GetYaxis().FindBin(cut)
#            eff_sig = histos['rel'+isoKey]['isovscen_rel'+tag].Integral(xbin+1,xbin+1,0,ybin)
#            eff_sig /= histos['rel'+isoKey]['isovscen_rel'+tag].Integral(xbin+1,xbin+1,0,nbinsy+1)
#            eff_bkg = histos['rel'+isoKey]['bkgisovscen_rel'+tag].Integral(xbin+1,xbin+1,0,ybin)
#            eff_bkg /= histos['rel'+isoKey]['bkgisovscen_rel'+tag].Integral(xbin+1,xbin+1,0,nbinsy+1)
#
#            effCenCurves[tag][0].SetPoint(xbin,cen,eff_sig)
#            effCenCurves[tag][1].SetPoint(xbin,cen,eff_bkg)
#            
#        
#    drawIsolationPlots(histos,channel)
#
#    return effCurves,effCenCurves
#
#
#
#
#
#def drawIsolationPlots(histos,channel):
#
#    """ draws the isolation histograms """
#   
#    c=ROOT.TCanvas('c','c',500,500)
#    c.SetTopMargin(0.05)
#    c.SetLeftMargin(0.12)
#    c.SetRightMargin(0.12)
#    c.SetBottomMargin(0.1)
#
#    #loop over isolation components
#    for iso in histos:
#
#        for key in ['rhovscen','isovsrho']:
#            if not key in histos[iso] : continue
#            c.Clear()
#            c.SetRightMargin(0.12)
#            histos[iso][key].Draw('colz')
#            histos[iso][key+'_prof'].SetMarkerStyle(20)
#            histos[iso][key+'_prof'].Draw('e1same')
#            canvasHeader(extraTxt=['#bf{[%s]}'%ISOTITLE[isoKey],histos[iso][key+'_prof_form']])
#            c.Modified()
#            c.Update()
#            c.RedrawAxis()
#            for ext in ['png','pdf']:
#                c.SaveAs('{0}_{1}.{2}'.format(channel,histos[iso][key].GetName(),ext))
#
#        for key in ['isovscen','isovscen_rel','isovscen_relccor'] :
#
#            if not key in histos[iso] : continue
#
#            c.SetRightMargin(0.03)
#            hsig,hbkg=histos[iso][key],histos[iso]['bkg'+key]
#            hsig.SetBarWidth(0.4);
#            hsig.SetBarOffset(-0.25);
#            hsig.SetLineWidth(2)
#            hbkg.SetBarWidth(0.4);
#            hbkg.SetBarOffset(0.25);
#            hbkg.SetLineColor(ROOT.kGreen+3);
#            hsig.Draw("candle1");
#            hbkg.Draw("candle1 same");            
#            leg=c.BuildLegend(0.73,0.94,0.93,0.8)
#            leg.SetBorderSize(0)
#            leg.SetTextFont(42)
#            leg.SetTextSize(0.035)
#            leg.SetFillStyle(0)
#
#            extraTxt=['#bf{[%s]}'%ISOTITLE[isoKey]]
#            if key+'_prof' in histos[iso]:
#                histos[iso][key+'_prof'].Draw('same')
#                extraTxt.append(histos[iso][key+'_prof_form'])
#            canvasHeader(extraTxt=extraTxt)
#            c.Modified()
#            c.Update()
#            c.RedrawAxis()
#            for ext in ['png','pdf']:
#                c.SaveAs('{0}_{1}.{2}'.format(channel,histos[iso][key].GetName(),ext))
#
#
#def plotIsoSummary( grColl, name, channel ):
#    
#    c=ROOT.TCanvas('c','c',500,500)
#    c.SetTopMargin(0.05)
#    c.SetLeftMargin(0.12)
#    c.SetRightMargin(0.03)
#    c.SetBottomMargin(0.1)
#    c.SetGridy()
#
#    mg=ROOT.TMultiGraph()
#    for g in grColl: mg.Add(g,'l')
#    mg.Draw('al')
#    if 'roc' in name:
#        mg.GetXaxis().SetTitle('Signal efficiency')
#        mg.GetYaxis().SetTitle('Background efficiency')
#        mg.GetXaxis().SetRangeUser(0,1)
#        mg.GetYaxis().SetRangeUser(0,1)
#        leg=c.BuildLegend(0.15,0.94,0.4,0.94-0.03*len(grColl))
#        l=ROOT.TLine(0.9,0,0.9,1)
#        l.SetLineColor(ROOT.kGray)
#        l.Draw()
#    elif 'effvscen' in name:
#        mg.GetXaxis().SetTitle('Centrality')
#        mg.GetYaxis().SetTitle('%s efficiency'%('Background' if 'bkg' in name else 'Signal'))
#        #mg.GetYaxis().SetRangeUser(0. if 'bkg' in name else 0.5,1)
#        #mg.GetXaxis().SetRangeUser(0,100)
#        leg=c.BuildLegend(0.65,0.7,0.92,0.7-0.03*len(grColl))
#    elif 'eff' in name:
#        mg.GetXaxis().SetTitle('max Relative Isolation')
#        mg.GetYaxis().SetTitle('%s efficiency'%('Background' if 'bkg' in name else 'Signal'))
#        mg.GetYaxis().SetRangeUser(0. if 'bkg' in name else 0.5,1)
#        mg.GetXaxis().SetRangeUser(-1,4)
#        leg=c.BuildLegend(0.65,0.7,0.92,0.7-0.03*len(grColl))
#    
#    leg.SetBorderSize(0)
#    leg.SetTextFont(42)
#    leg.SetTextSize(0.035)
#    leg.SetFillStyle(0)   
#    canvasHeader()
#    c.Modified()
#    c.Update()
#    c.RedrawAxis()
#    for ext in ['png','pdf']:
#        c.SaveAs("%s_isosummary%s.png"%(channel,name))
#
#
##open all data files which are relevant
#url=sys.argv[1]
#channel=sys.argv[2]
#dataTag='SkimElectrons' if channel=='ee' else 'SkimMuons'
#dataCat=121 if channel=='ee' else 169
#t=ROOT.TChain('tree')
#for f in os.listdir(url):
#    if not dataTag in f or not '.root' in f: continue
#    t.Add(os.path.join(url,f))
#
##loop over events and filter the required channel
#leptons=[]
#for iev in range(t.GetEntries()):
#    t.GetEntry(iev)
#
#    #get dilepton
#    lepColl=getLeptons(t)
#    dil=dileptonBuilder(lepColl)
#    if not dil.isZ and not dil.isSS : continue 
#    if (dil.l1.pdgId*dil.l1.pdgId)!=dataCat : continue
#    if dil.p4.M()<20 : continue
#
#    #centrality
#    cenbin=t.cenbin
#    leptons.append( [dil.l1,dil.isZ,dil.isSS,cenbin] )
#    leptons.append( [dil.l2,dil.isZ,dil.isSS,cenbin] )
#
#print 'Found',len(leptons),'leptons for optimization'
#
##centrality bins
#zCentralities=[]
#for _,isZ,_,cenbin in leptons:
#    if not isZ: continue
#    zCentralities.append(cenbin)
#zcenq=np.percentile(zCentralities,[0,20,40,60,80,100],axis=0)
#
#
#rocCurves=[]
#effVsCen=[]
#for isoKey in ['iso_rhosub','iso20_rhosub','iso25_rhosub','iso30_rhosub','miniiso_rhosub']: #,'iso_nosub','iso_rhosub','iso_partrhosub','iso_leprhosub','chiso_nosub','chiso_chrhosub']:
#    irocs,ieffvscen=tuneIsolation(isoKey=isoKey,isoComponents=ISOCOLLECTION[isoKey],leptons=leptons,channel=channel)
#    rocCurves.append( (isoKey,        irocs['']) )
#    #rocCurves.append( (isoKey+'ccor', irocs['ccor']) )
#
#    effVsCen.append( (isoKey,         ieffvscen['']) )
#    #effVsCen.append( (isoKey+'ccor',  ieffvscen['ccor']) )
#
#
#
#plotIsoSummary( grColl=[x[0] for _,x in rocCurves], name='roc',    channel=channel )
#plotIsoSummary( grColl=[x[1] for _,x in rocCurves], name='sigeff', channel=channel )
#plotIsoSummary( grColl=[x[2] for _,x in rocCurves], name='bkgeff', channel=channel )
#
#plotIsoSummary( grColl=[x[0] for _,x in effVsCen], name='sigeffvscen', channel=channel )
#plotIsoSummary( grColl=[x[1] for _,x in effVsCen], name='bkgeffvscen', channel=channel )
#
#print 'Working points'
#print 'Isolation\t Efficiency\t Cut'
#for i in range(len(rocCurves)):
#    print '%s\t%3.3f\t%3.3f'%(rocCurves[i][0],rocCurves[i][1][3],rocCurves[i][1][4])
