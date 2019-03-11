import ROOT
import os,sys
import numpy as np
from HeavyIonsAnalysis.topskim.Plot import *

LUMI=1618.466*(1e-6)
ABEAM=208
SAMPLES=[
    ('WJetsToLNu_TuneCP5_5020GeV-amcatnloFXFX-pythia8',        'W',                    21159*(ABEAM**2)*LUMI,          17),         
    ('DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8', 'Z/#gamma^{*}',         2010*(ABEAM**2)*LUMI,           "#fdc086"),         
    ('Skim',                                                   'Z/#gamma^{*} (data)',  None,                           "#fdc086"),
    ('TT_TuneCP5_5p02TeV-powheg-pythia8',                      't#bar{t}',             69*(ABEAM**2)*LUMI,             633),
    ('ST_tW_antitop_5f_NoFullyHadronicDecays_hdampDOWN_TuneCP5_5p02TeV-powheg-pythia8', 'tW', 3.04*(ABEAM**2)*LUMI,    "#7fc97f"),
    ('ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5down_5p02TeV-powheg-pythia8',           'tW', 3.04*(ABEAM**2)*LUMI,    "#7fc97f"),
    ('WWTo2L2Nu_NNPDF31_TuneCP5_5p02TeV-powheg-pythia8',       'VV',                   1.21*(ABEAM**2)*LUMI,           "#386cb0"),         
    ('WZTo3LNU_NNPDF30_TuneCP5_5p20TeV-powheg',                'VV',                   1.77*(ABEAM**2)*LUMI,           "#386cb0"),         
    ('Skim',                                                   'Combinatorial (data)', None,                           17),
    ('Skim',                                                   'Data',                 None,                           1),
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
                print e
                pass
        fIn.Close()

    return histos

def showRateVsRun(url,catList=['e','m']):

    """ checks for outliers in observed rate """

    histos=getDataSummedUp(url)
    p=Plot('ratevsrun',com='5.02 TeV')
    p.spimposeWithErrors=True
    p.forceRange=[0,5000]
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


def makeControlPlot(url,cat,pname,dyFromData=True,combFromData=True):

    """makes a standard control plot using the ROOT files in 'url' 
    'cat' and 'pname' are used to build the full plot name
    'dyFromData' and 'combFromData' are used to set flag 
    if these processes are taken from the data-driven estimates """

    plotsPerProc={}
    for tag,title,scale,color in SAMPLES:

        icat=cat
        if title=='Combinatorial (data)': icat='ss'+cat
        if title=='Z/#gamma^{*} (data)': continue
        normByWgtSum=True if scale else False

        plots=getDataSummedUp(url,[icat],pname,tag,normByWgtSum)
        if not icat in plots: continue
        h=plots[icat]
        if scale : h.Scale(scale)

        if not title in plotsPerProc:
            plotsPerProc[title]=h.Clone('%s_%s'%(pname,title))               
            plotsPerProc[title].SetTitle(title)
            ci=ROOT.TColor.GetColor(color) if isinstance(color,str) else color
            plotsPerProc[title].SetLineColor(ci)
            plotsPerProc[title].Reset('ICE')
        plotsPerProc[title].Add(h)

    p=Plot('%s_%s'%(cat,pname) ,com='5.02 TeV')

    for proc in ['Z/#gamma^{*}', 
                 'Combinatorial (data)' if combFromData else 'W',                 
                 't#bar{t}',
                 'tW', 
                 'VV', 
                 'Data']:
        if not proc in plotsPerProc : continue        
        p.add(h=plotsPerProc[proc],
              title=proc,
              color=plotsPerProc[proc].GetLineColor(),
              isData=True if proc=='Data' else False,
              spImpose=False,
              isSyst=False)
    p.savelog=True
    p.show(outDir='./',lumi=LUMI)
    p.reset()


url=sys.argv[1]

showRateVsRun(url)

cats=[]
cats+=['zee','zmm']
cats+=['zeeBB','zeeEE']
#cats+=['zeeZawaypfj']
for cat in cats:
    for d in ['mll','ptll','dphill', 'detall',
              'l1pt','l1eta','l1chreliso','l1phoreliso','l1neureliso',
              'l2pt','l2eta','l2chreliso','l2phoreliso','l2neureliso',
              'npfjets','npfbjets','pf1jpt','pf1jeta','pf1jcsv','pf2jpt','pf2jeta','pf2jcsv',
              #'ntkjets','ntkbjets','tk1jbalance','tk1jpt','tk1jeta','tk1jcsv','tk2jbalance','tk2jpt','tk2jeta','tk2jcsv','tkrho'
              ]:
        makeControlPlot(url,cat,d,dyFromData=True,combFromData=True)
