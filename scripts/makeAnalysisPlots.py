import ROOT
import os,sys
import numpy as np
from HeavyIonsAnalysis.topskim.Plot import *
from HeavyIonsAnalysis.topskim.RinRout import *

LUMI=1618.466*(1e-6)
CHLUMI={'mm':1587.941*(1e-6),
        'ee':1664.148*(1e-6),
        'blind':446.931*(1e-6)}

ABEAM=208
SAMPLES=[
    ('WJetsToLNu_TuneCP5_5020GeV-amcatnloFXFX-pythia8',          'W',                    21159*(ABEAM**2)*LUMI,        17),         
    ('DYJetsToLL_M-10to50_TuneCP5_5020GeV-amcatnloFXFX-pythia8', 'Z/#gamma^{*}',         1506*(ABEAM**2)*LUMI,         "#fdc086"),         
    ('DYJetsToLL_TuneCUETP8M1_5020GeV-amcatnloFXFX-pythia8',     'Z/#gamma^{*}',         2010*(ABEAM**2)*LUMI,         "#fdc086"),         
    ('TT_TuneCP5_5p02TeV-powheg-pythia8',                        't#bar{t}',             69*(ABEAM**2)*LUMI,           633),
    ('ST_tW_antitop_5f_NoFullyHadronicDecays_hdampDOWN_TuneCP5_5p02TeV-powheg-pythia8', 'tW', 3.04*(ABEAM**2)*LUMI,    "#7fc97f"),
    ('ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5down_5p02TeV-powheg-pythia8',           'tW', 3.04*(ABEAM**2)*LUMI,    "#7fc97f"),
    ('WWTo2L2Nu_NNPDF31_TuneCP5_5p02TeV-powheg-pythia8',         'VV',                   1.21*(ABEAM**2)*LUMI,         "#386cb0"),         
    ('WZTo3LNU_NNPDF30_TuneCP5_5p20TeV-powheg',                  'VV',                   1.77*(ABEAM**2)*LUMI,         "#386cb0"),         
    ('ZZTo4L_5p02TeV_powheg_pythia8',                            'VV',                   0.441*(ABEAM**2)*LUMI,        "#386cb0"),         
    ('Skim',                                                     'Comb. (data)',         None,                         17),
    ('Skim',                                                     'Data',                 None,                         1),
]


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)


def getDataSummedUp(url,catList=['e','m'],pname='ratevsrun',tag='Skim',normByWgtSum=False):

    """
    loops over a directory and sums all plots 'pname' for different categories in 'catList' 
    'tag' is used to match the start of the files name, for the plots to be added
    'normByWgtSum' is used to read the weight sum and scale the histograms to return
    """


    #sum up plots from data
    histos={}
    for f in os.listdir(url):

        #open file
        if f.find(tag)!=0 : continue
        fIn=ROOT.TFile.Open(os.path.join(url,f))

        scale=1.0
        if normByWgtSum:
            scale=fIn.Get('wgtsum').GetBinContent(1)

        #read histogram
        for c in catList:
            h=fIn.Get('%s_%s'%(c,pname))        
            try:
                h.Scale(1./scale)
                if not c in histos: 
                    histos[c]=h.Clone()
                    histos[c].SetDirectory(0)
                else:
                    histos[c].Add(h)
            except Exception as e:
                #print e
                pass
        fIn.Close()

    return histos

def showRateVsRun(url,catList=['e','m']):

    """ checks for outliers in observed rate """

    histos=getDataSummedUp(url)
    p=Plot('ratevsrun',com='5.02 TeV')
    p.spimposeWithErrors=True
    p.forceRange=[0,5000]
    p.plotformats = ['pdf','png','root']
    for c in catList:

        rates=np.array([[i+1,histos[c].GetBinContent(i+1)] for i in range(histos[c].GetNbinsX())])
        q=np.percentile(rates,50,axis=0)[1]
        lowRates=rates[rates[:,1]<q*0.5][:,0].flatten()
        runs=[int(histos[c].GetXaxis().GetBinLabel(int(i))) for i in lowRates]        
        print 'Runs for',c,'with rate < %f/2'%q
        print runs

        #add to plot
        title='e' if c[0]=='e' else '#mu'
        color=2 if c[0]=='e' else 8
        marker=20 if c[0]=='e' else 24
        histos[c].SetMarkerStyle(marker)
        histos[c].SetMarkerColor(color)
        histos[c].SetTitle(title)        
        p.add(histos[c],title=title,color=color,isData=False,spImpose=True,isSyst=False)
    p.show(outDir='./',lumi=LUMI)


def makeControlPlot(url,cat,pname,dyFromData,combFromData,dySF,plotter=None,rebin=None,saveTeX=False):

    """makes a standard control plot using the ROOT files in 'url' 
    'cat' and 'pname' are used to build the full plot name
    'dyFromData' and 'combFromData' are used to set flag 
    if these processes are taken from the data-driven estimates """

    plotsPerProc={}
    for tag,title,scale,color in SAMPLES:

        icat=cat
        ilumi=LUMI
        if 'mm' in icat : ilumi=CHLUMI['mm']
        if 'ee' in icat : ilumi=CHLUMI['ee']
        if not 'z' in icat: ilumi=CHLUMI['blind']
        if title=='Comb. (data)': icat='ss'+cat
        normByWgtSum=True if scale else False

        plots=getDataSummedUp(url,[icat],pname,tag,normByWgtSum)
        if not icat in plots: continue
        h=plots[icat]
        if scale :

            scale *= ilumi/LUMI

            #special treatment for DY (scale and shape)
            zShapeCat=''
            if dyFromData and title=='Z/#gamma^{*}':                                              
                if icat.find('mm') in [0,1]: 
                    scale*=dySF['mm'][0]
                    zShapeCat='zmm' if icat[0]=='m' else ''
                if icat.find('ee') in [0,1]: 
                    scale*=dySF['ee'][0]
                    zShapeCat='zee' if icat[0]=='e' else ''
                if icat.find('em') in [0,1]: 
                    scale*=dySF['em'][0]
                    zShapeCat='zmm' if icat[0]=='e' else ''

            h.Scale(scale)
            
            #change shape
            if dyFromData>1 and len(zShapeCat)>0:
                zplot=getDataSummedUp(url,[zShapeCat],pname,'Skim',False)[zShapeCat]
                sszplot=getDataSummedUp(url,['ss'+zShapeCat],pname,'Skim',False)['ss'+zShapeCat]
                zplot.Add(sszplot,-1)
                totz=zplot.Integral()
                if totz>0:
                    zplot.Scale(h.Integral()/totz)
                    h.Reset('ICE')
                    h.Add(zplot)
                zplot.Delete()
                sszplot.Delete()

        plots=getDataSummedUp(url,[icat],pname,tag,normByWgtSum)

        if title=='Comb. (data)' and not 'z' in icat:
            h.Scale(CHLUMI['blind']/LUMI)

        if not title in plotsPerProc:
            plotsPerProc[title]=h.Clone('%s_%s'%(pname,title))               
            plotsPerProc[title].SetTitle(title)
            ci=ROOT.TColor.GetColor(color) if isinstance(color,str) else color
            plotsPerProc[title].SetLineColor(ci)
            plotsPerProc[title].Reset('ICE')
        plotsPerProc[title].Add(h)

    p=Plot('%s_%s'%(cat,pname) ,com='5.02 TeV')

    for proc in ['Z/#gamma^{*}', 
                 'Comb. (data)' if combFromData else 'W',                 
                 't#bar{t}',
                 'tW', 
                 'VV', 
                 'Data']:
        if not proc in plotsPerProc : continue        
        
        if rebin: plotsPerProc[proc].Rebin(rebin)
        
        p.add(h=plotsPerProc[proc],
              title=proc,
              color=plotsPerProc[proc].GetLineColor(),
              isData=True if proc=='Data' else False,
              spImpose=False,
              isSyst=False)
    p.savelog=True
    p.show(outDir='./',lumi=ilumi,saveTeX=saveTeX if saveTeX else False)
    if plotter : p.appendTo(plotter)
    p.reset()

def computeDYScaleFactors(url):

    cats=['zee','zmm','ee','mm','em']
    data=getDataSummedUp(url,cats,'mll','Skim',False)
    dy=getDataSummedUp(url,cats,'mll','DYJetsToLL_TuneCUETP8M1_5020GeV-amcatnloFXFX-pythia8',True)
    nin={}
    nout={}
    ndyin={}
    ndyout={}
    intErr=ROOT.Double(0)
    for c in ['ee','mm','em']:        
        if c != 'em':
            data[c].Add(data['z'+c])
            dy[c].Add(dy['z'+c])
            dy[c].Scale(2010*(ABEAM**2)*LUMI)
            
        binMin,binMax=data[c].GetXaxis().FindBin(20),data[c].GetNbinsX()+1
        bin1,bin2=data[c].GetXaxis().FindBin(76),data[c].GetXaxis().FindBin(106)

        nin[c]=(data[c].Integral(bin1,bin2),0.)
        nout[c]=(data[c].Integral(binMin,binMax),0.)

        dyCts=dy[c].IntegralAndError(bin1,bin2,intErr)
        ndyin[c]=(dyCts,intErr)
        dyCts=dy[c].IntegralAndError(binMin,binMax,intErr)
        ndyout[c]=(dyCts,intErr)


    #table a la TOP-16-015 
    dySF={}
    dySF['ee']=doRinRout(nin,nout,ndyin,ndyout,'ee','mm')
    dySF['mm']=doRinRout(nin,nout,ndyin,ndyout,'mm','ee')
    sfProd=dySF['ee'][0]*dySF['mm'][0]
    sfProdUnc=ROOT.TMath.Sqrt((dySF['ee'][0]*dySF['mm'][1])**2+(dySF['mm'][0]*dySF['ee'][1])**2)        
    dySF['em']=(ROOT.TMath.Sqrt(sfProd),0.5*sfProdUnc/ROOT.TMath.Sqrt(sfProd))
    print 'SF_{em}=%f +/- %f'%(dySF['em'][0],dySF['em'][1])
    return dySF

def compareElectrons(url,dist):

    cats=['zee','zeeBB','zeeEB','zeeEE']
    cats+=['zeeafter','zeeBBafter','zeeEBafter','zeeEEafter']
    cats+=['zeebefore','zeeBBbefore','zeeEBbefore','zeeEEbefore']    
    data=getDataSummedUp(url,cats,dist,'Skim',False)
    for key in data: 
        tot=data[key].Integral()
        data[key].Scale(1./tot)
        data[key].GetYaxis().SetTitle('PDF')

    for pname,triplet,title in [('period',['zee','zeeafter','zeebefore'],['inclusive','run#geq327402','run<327402']),
                                ('region',['zee','zeeBB','zeeEB','zeeEE'],['inclusive','barrel-barrel','barrel-endcap','endcap-endcap'])]:

        p=Plot('%s_%s'%(dist,pname),com='5.02 TeV')
        p.doPoissonErrorBars=False
        p.add(data[triplet[0]],title=title[0],color=1,isData=True,spImpose=False,isSyst=False)
        p.add(data[triplet[1]],title=title[1],color=2,isData=False,spImpose=False,isSyst=False)
        p.add(data[triplet[2]],title=title[2],color=8,isData=False,spImpose=False,isSyst=False)
        if len(triplet)>3:
            p.add(data[triplet[3]],title=title[3],color=6,isData=False,spImpose=False,isSyst=False)
        p.show(outDir='./',lumi=LUMI,noStack=True)

def doEleIDPlots(url):

    for reg in ['EB','EE']:
        cats=['zeectrl'+reg,'sszeectrl'+reg]
        for dist in ['esihih','edetavtx','edphivtx','ehoe','eempinv','ed0', 'edz']:
            plots=getDataSummedUp(url,cats,dist,'Skim',False)
            
            #show the plots for simple variables
            p=Plot('%s%s'%(dist,reg),com='5.02 TeV')
            p.add(plots[cats[0]],title='Z#rightarrowee (data)',color=1,isData=True,spImpose=False,isSyst=False)
            p.add(plots[cats[1]],title='Comb. (data)',color=17,isData=False,spImpose=False,isSyst=False)
            p.savelog=True
            p.show(outDir='./',lumi=LUMI)

def doMuIDPlots(url):

    cats=['zmmctrl','sszmmctrl']
    for dist in ['mmusta', 'mtrklay', 'mchi2ndf', 'mmuhits', 'mpxhits', 'md0','mdz']:
        plots=getDataSummedUp(url,cats,dist,'Skim',False)
            
        #show the plots for simple variables
        p=Plot('%s'%dist,com='5.02 TeV')
        p.add(plots[cats[0]],title='Z#rightarrow#mu#mu (data)',color=1,isData=True,spImpose=False,isSyst=False)
        p.add(plots[cats[1]],title='Comb. (data)',color=17,isData=False,spImpose=False,isSyst=False)
        p.savelog=True
        p.show(outDir='./',lumi=LUMI)
        

def doIsolationROCs(url,ch='ee'):

    cats=['z'+ch,'ssz'+ch,'ss'+ch]

    #read all isolation plots and sum the contributions of the two leptons
    data={}
    for c in ['ch','pho','nh','']:
        for d in ['isovscen','isovsrho','isopvscen','iso','reliso','isop','relisop']:
            dist=c+d 
            for l in ['l1','l2']:
                plots=getDataSummedUp(url,cats,l+dist,'Skim',False)
                if l=='l1':
                    data[dist]=plots
                else:
                    for key in data[dist]:                           
                        data[dist][key].Add( plots[key] )

    cnv=ROOT.TCanvas('c','c',500,500)
    cnv.SetLeftMargin(0.12)
    cnv.SetTopMargin(0.05)
    cnv.SetBottomMargin(0.11)
    cnv.SetRightMargin(0.12)
    for c in ['ch','pho','nh','']:
        for d in ['isovscen','isovsrho','isopvscen']:
            iso=c+d
            #subtract combinatorial background
            if not iso in data: continue
            if not 'z'+ch in data[iso]: continue
            sig=data[iso]['z'+ch].Clone('sig'+iso)
            sig.Add(data[iso]['ssz'+ch],-1)
            px=sig.ProfileX()
            px.SetMarkerStyle(20)            
            if not 'isop' in iso and 'vsrho' in iso:
                func=ROOT.TF1('func','[0]*([1]/([2]+x)+1.0)',0,2)                
                px.Fit(func)
            sig.Draw('colz')
            px.Draw('e1same')
            txt=ROOT.TLatex()
            txt.SetNDC(True)
            txt.SetTextFont(42)
            txt.SetTextSize(0.045)
            txt.SetTextAlign(12)
            txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
            if not 'isop' in iso and 'vsrho' in iso:
                txt.DrawLatex(0.15,0.88,
                              '<UE>(#rho)=%3.3f#times(1+#frac{%3.3f}{%3.3f+#rho})'%(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2)))                
            cnv.Modified()
            cnv.Update()
            for ext in ['png','pdf']:
                cnv.SaveAs('%s.%s'%(iso,ext))

                
    #add the two leptons and subtract to the signal
    rocs=[]
                
    isoDists={
        'iso':(1,'I'),
        'reliso':(4,"I_{rel}"),
        'isop':(2,"I'"),
        'relisop':(8,"I_{rel}'"),
        }
    for iso in ['iso','reliso','isop','relisop']:

        sig=data[iso]['z'+ch].Clone('sig'+iso)
        sig.Add(data[iso]['ssz'+ch],-1)
        bkg=data[iso]['ss'+ch].Clone('bkg'+iso)
        sig.Rebin()
        bkg.Rebin()
        color,title=isoDists[iso]
    
        rocs.append( ROOT.TGraph() )
        tot_sig=sig.Integral()
        tot_bkg=bkg.Integral()
        nbins=sig.GetNbinsX()
        best_effsig=0
        best_xbin=nbins+1
        for xbin in range(nbins+1):
            effsig=sig.Integral(0,xbin+1)/tot_sig
            rocs[-1].SetPoint(xbin,effsig,bkg.Integral(0,xbin+1)/tot_bkg)
            if abs(effsig-0.9)>abs(best_effsig-0.9) : continue
            best_effsig=effsig
            best_xbin=xbin
        print iso,best_effsig,best_xbin,sig.GetXaxis().GetBinCenter(best_xbin+1)
        rocs[-1].SetTitle(title)
        rocs[-1].SetLineColor(color)
        rocs[-1].SetMarkerColor(color)
        rocs[-1].SetLineWidth(2)

    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.11)
    mg=ROOT.TMultiGraph()
    for g in rocs: mg.Add(g,'l')
    mg.Draw('al')
    mg.GetXaxis().SetRangeUser(0,1)
    mg.GetYaxis().SetRangeUser(0,1)
    mg.GetXaxis().SetTitle('Signal efficency')
    mg.GetYaxis().SetTitle('Background efficency')
    leg=c.BuildLegend(0.15,0.94,0.4,0.6)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.05)
    leg.SetFillStyle(0)
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    for ext in ['png','pdf']:
        c.SaveAs('isorocs_%s.%s'%(ch,ext))



def doJetHotSpots(url,cats):

    data={}
    for dist in ['pf1jetavsphi','pf2jetavsphi']:
        data[dist]=getDataSummedUp(url,cats,dist,'Skim',False)

    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.11)
    c.SetRightMargin(0.12)
    h=data['pf1jetavsphi'][cats[0]].Clone('etavsphi')
    h.Reset('ICE')
    for cat in cats:
        h.Add(data['pf1jetavsphi'][cat])
        h.Add(data['pf2jetavsphi'][cat])
        print cat,h.Integral()
    h.Draw('colz')
    h.Rebin2D(2,2)
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    for ext in ['png','pdf']:
        c.SaveAs('jetetavsphi.%s'%ext)


def checkAcceptance(url):

    cats=['zee','zeehpur','zeehpurBB']
    zee=getDataSummedUp(url,cats,'mll','Skim',False)
    sszee=getDataSummedUp(url,['ss'+x for x in cats],'mll','Skim',False)
    ttcats=[x[1:] for x in cats]
    ttbar=getDataSummedUp(url,ttcats,'mll','TT_TuneCP5_5p02TeV-powheg-pythia8',True)
    nzee_ini=0
    nsszee_ini=0
    ntt_ini=0
    for i in range(len(cats)):
        nzee=zee[cats[i]].Integral()
        nsszee=sszee['ss'+cats[i]].Integral()
        nzee-=nsszee
        ntt=ttbar[ttcats[i]].Integral()
        if i==0: 
            nzee_ini=nzee
            nsszee_ini=nsszee
            ntt_ini=ntt
        print cats[i],nzee/nzee_ini,nsszee/nsszee_ini,ntt/ntt_ini


url=sys.argv[1]
#doEleIDPlots(url)
#doIsolationROCs(url,'ee')
#doMuIDPlots(url)
#doIsolationROCs(url,'mm')

#showRateVsRun(url)

dySF=computeDYScaleFactors(url)

#compareElectrons(url,'mll')
#compareElectrons(url,'ptll')
#compareElectrons(url,'detall')



#doJetHotSpots(url,['zmm','zee','em'])

cats=[]
cats+=['zee','zmm','mm','em','ee']
cats+=['zeehpur','zmmhpur']
cats+=['mm','em','ee']
cats+=['mmhpur','emhpur','eehpur']
cats+=['mm0pfb','mmgeq1pfb','em0pfb','emgeq1pfb','ee0pfb','eegeq1pfb',]
cats+=['mmhpur0pfb','mmhpurgeq1pfb','emhpur0pfb','emhpurgeq1pfb','eehpur0pfb','eehpurgeq1pfb',]
for cat in cats:
    for d in ['mll','ptll','l1pt','l1eta','l2pt','l2eta']:                        
        makeControlPlot(url,cat,d,1,True,dySF)

fIn=ROOT.TFile.Open('plotter.root','RECREATE')
fIn.Close()
for cat in cats:
    for d in ['acopl', 'detall','drll']:
        makeControlPlot(url,cat,d,1,True,dySF,'plotter.root',rebin=2) #,saveTeX=True)


for cat in cats:
    for d in ['npfjets','npfbjets','pf1jpt','pf1jeta','pf1jcsv','pf2jpt','pf2jeta','pf2jcsv','pfrapavg','pfraprms','pfrapmaxspan']:
        makeControlPlot(url,cat,d,2,True,dySF)

for cat in cats:
    for d in ['pfht','pfmht']:        
        makeControlPlot(url,cat,d,2,True,dySF,rebin=2)
              

#checkAcceptance(url)
