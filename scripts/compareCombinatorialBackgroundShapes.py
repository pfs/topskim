from HeavyIonsAnalysis.topskim.Plot import *

LUMI=1618.466*(1e-6)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def getDistribution(pname='mll',flav='ee',pfix='',url='combbackground_plots.root'):

    """
    compares the distributions in a given flavour
    """
    cats=[ (flav              , 'SS',                 ROOT.kBlack),            
           (flav+'iso'        , 'iso SS',             ROOT.kRed),
           (flav+'mix'        , 'mix',                ROOT.kMagenta),
           (flav+'mixwgt'     , 'mix (wgt)',          ROOT.kAzure+1),
           (flav+'mixiso'     , 'iso mix',            ROOT.kYellow-3),
           (flav+'mixwgtiso'  , 'iso mix (wgt)',      ROOT.kYellow+3)
           ]

    fIn=ROOT.TFile.Open(url)
    histos=[]
    histosNorm={}
    for cat,title,ci in cats:
        key=pname+'_'+cat+pfix
        histos.append( fIn.Get(key) )        
        try:
            histos[-1].SetDirectory(0)
            histos[-1].SetTitle(title)
            histos[-1].SetLineColor(ci)
            histosNorm[cat]=histos[-1].Integral()
            if not 'mix' in cat: continue
            normCat=cat.replace('mix','')
            normCat=normCat.replace('wgt','')
            sf = histosNorm[normCat]/histosNorm[cat]
            histos[-1].Scale(sf)
        except Exception as e:
            histos.pop(1)
    return histos

url='comb_shapes_plotter.root'
fIn=ROOT.TFile.Open(url,'RECREATE')
fIn.Close()
for flav in ['ee','em','mm'] :
    for pname in ['mll','ptll','acopl','lpt','leta','liso'] :
        pfixes=['']
        if pname[0]=='l' : pfixes=['lead','sublead']
        for pfix in pfixes:
            histos=getDistribution(pname,flav,pfix)
            p=Plot('%s_%s%s'%(flav,pname,pfix),com='#sqrt{s_{NN}}=5.02 TeV')
            #p.spimposeWithErrors=True
            for h in histos:
                p.add(h,title=h.GetTitle(),color=h.GetLineColor(),isData=False,spImpose=True,isSyst=False)
            p.show(outDir='./',lumi=LUMI)
            p.appendTo(url)
