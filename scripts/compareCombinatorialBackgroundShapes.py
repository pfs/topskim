from HeavyIonsAnalysis.topskim.Plot import *

LUMI=1618.466*(1e-6)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def getDistribution(pname='mll',flav='ee',pfix='',url='combbackground_plots_data.root',isMC=False):

    """
    compares the distributions in a given flavour
    """
    cats=[ (flav              , 'SS',                 ROOT.kBlack),            
           #(flav+'iso'        , 'iso SS',             ROOT.kRed),
           (flav+'mix'        , 'mix',                ROOT.kMagenta),
           #(flav+'mixiso'     , 'iso mix',            ROOT.kYellow-3),
           ]
    if isMC:
        cats=[
            (flav              , 'OS',                 ROOT.kGray),            
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
            if len(histos):
                histos.pop(len(histos)-1)
    return histos

url='comb_shapes_plotter.root'
fIn=ROOT.TFile.Open(url,'RECREATE')
fIn.Close()
for flav in ['ee','em','mm'] :
    for pname in ['mll','ptll','acopl','sphericity','lpt','leta','liso'] :
        pfixes=['']
        if pname[0]=='l' : pfixes=['lead','sublead']
        for pfix in pfixes:
            histos=getDistribution(pname,flav,pfix)
            histos_w=getDistribution(pname,flav,pfix,url='combbackground_plots_w.root',isMC=True)

            p=Plot('%s_%s%s'%(flav,pname,pfix),com='#sqrt{s_{NN}}=5.02 TeV')
            #p.spimposeWithErrors=True
            for h in histos:
                p.add(h,title=h.GetTitle(),color=h.GetLineColor(),isData=False,spImpose=True,isSyst=False)
            try:
                for h in histos_w:
                    h.Scale( histos[0].Integral()/h.Integral())
                    p.add(h,title=h.GetTitle()+' (W)',color=h.GetLineColor(),isData=False,spImpose=True,isSyst=False)
            except:
                pass
            p.show(outDir='./',lumi=LUMI)
            p.appendTo(url)
