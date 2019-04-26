#!/usr/bin/env python

import ROOT
import sys
import os
import pickle
import random
from array import array
from collections import defaultdict
from HeavyIonsAnalysis.topskim.EventReader import *

BDTMETHOD=ROOT.TString('BDTG')
BDTWGTS="/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/HeavyIonsAnalysis/topskim/scripts/training_dy/weights/TMVAClassification_BDTG.weights.xml"

def prepareDileptonCollection(url):

    """loops over all the available events and stores the information on the dileptons in each event"""

    #build the chain
    t=ROOT.TChain('tree')
    for f in os.listdir(url):
        if not '.root' in f: continue
        if not 'Skim' in f : continue
        t.Add(os.path.join(url,f))

    #loop over events
    print 'Analysing',t.GetEntries(),'events'
    dilCollection=defaultdict(list)
    for iev in range(t.GetEntries()):
        t.GetEntry(iev)
        try:
            dil=getDilepton(t,[13,11])
            dilCollection[(dil.flavour,dil.isZ)].append(dil)
        except Exception as e:
            print e
            pass

    #save dict in a cache
    pckURL='dilepton_summary.pck'
    with open(pckURL,'w') as cache:
        pickle.dump( dilCollection,cache,pickle.HIGHEST_PROTOCOL)
    return pckURL

def createMixedFriendTrees(url,mixFile,outURL):
    
    """for each file it mixes one lepton at each time and dumps a friend tree with mix_lep_* branches"""

    os.system('mkdir -p %s'%outURL)

    #read the dilepton collection to use for the mixing
    with open(mixFile,'r') as cache:
        mixDileptons=pickle.load(cache)

    #start tmva reader
    tmva_reader = ROOT.TMVA.Reader( "!Color:!Silent" );
    bdt_l1pt=array('f',[0])
    tmva_reader.AddVariable('lep_pt[0]',bdt_l1pt)
    bdt_apt=array('f',[0])
    tmva_reader.AddVariable('apt',bdt_apt)
    bdt_llpt=array('f',[0])
    tmva_reader.AddVariable('llpt',bdt_llpt)
    bdt_abslleta=array('f',[0])
    tmva_reader.AddVariable('abs(lleta)',bdt_abslleta)
    bdt_dphill=array('f',[0])
    tmva_reader.AddVariable('abs(dphi)',bdt_dphill)
    bdt_sumabseta=array('f',[0])
    tmva_reader.AddVariable('abs(lep_eta[0])+abs(lep_eta[1])',bdt_sumabseta)
    tmva_reader.BookMVA(BDTMETHOD,BDTWGTS)

    #mix each file separately
    for f in os.listdir(url):
        if not '.root' in f: continue
        if not 'Skim' in f : continue

        out_f=ROOT.TFile.Open(os.path.join(outURL,f),'RECREATE')
        out_t=ROOT.TTree('tree','tree')
        out_t_branches={}
        for name in LEPTONBRANCHES:
            out_t_branches[name]=array('f',[0.,0.])
            out_t.Branch('lep_'+name,out_t_branches[name],'lep_%s[2]/F'%name)
        for name in DILEPTONBRANCHES:
            out_t_branches[name]=array('f',[0.])
            out_t.Branch(name,out_t_branches[name],'%s/F'%name)

        orig_t=ROOT.TChain('tree')
        orig_t.Add(os.path.join(url,f))

        print 'Mixing',orig_t.GetEntries(),'events from',f,'output @',out_f.GetName()
        for iev in range(orig_t.GetEntries()):

            #original event
            orig_t.GetEntry(iev)
            orig_dil  = getDilepton(orig_t,[13,11])
            orig_flav = orig_dil.flavour
            orig_isZ  = orig_dil.isZ

            #mix one lepton at each time
            for i in range(2):
                leptons=[getattr(orig_dil,'l%d'%(i+1))]
                lidToGet=orig_dil.l2.pdgId if i==0 else orig_dil.l1.pdgId
                while True:
                    mix_dil=random.choice(mixDileptons[orig_flav,orig_isZ])
                    leptons.append( mix_dil.l1 if mix_dil.l1.pdgId==lidToGet else mix_dil.l2 )
                    break

                #sort by pT
                leptons.sort(key=lambda x: x.pt, reverse=True)

                #fill the out tree
                for name in LEPTONBRANCHES:
                    for il in range(2):                        
                        out_t_branches[name][il]=getattr(leptons[il],name)
                
                llp4             = leptons[0].p4+leptons[1].p4
                dphill           = abs(leptons[0].p4.DeltaPhi(leptons[1].p4))
                detall           = abs(leptons[0].eta-leptons[1].eta)
                sumeta           = leptons[0].eta+leptons[1].eta
                bdt_l1pt[0]      = max(leptons[0].pt,leptons[1].pt)
                bdt_llpt[0]      = llp4.Pt()
                bdt_apt[0]       = abs(leptons[0].pt-leptons[1].pt)/(leptons[0].pt+leptons[1].pt)
                bdt_abslleta[0]  = abs(llp4.Eta())
                bdt_dphill[0]    = dphill
                bdt_sumabseta[0] = abs(leptons[0].eta)+abs(leptons[1].eta)                
                out_t_branches['llpt'][0]       = llp4.Pt()
                out_t_branches['lleta'][0]      = llp4.Eta()
                out_t_branches['llphi'][0]      = llp4.Phi()
                out_t_branches['llm'][0]        = llp4.M()
                out_t_branches['dphi'][0]       = dphill
                out_t_branches['deta'][0]       = detall
                out_t_branches['sumeta'][0]     = sumeta
                out_t_branches['apt'][0]        = bdt_apt[0]
                out_t_branches['bdt'][0]        = tmva_reader.EvaluateMVA(BDTMETHOD)
                out_t_branches['bdtrarity'][0]  = tmva_reader.GetRarity(BDTMETHOD)

                out_t.Fill()
                    
            
        #write tree
        out_f.cd()
        out_t.SetDirectory(out_f)
        out_t.Write()
        out_f.Close()


def main():

    if len(sys.argv)<3:
        print 'python prepareCombinatorialBackgroundTree.py input_directory output_directory'
        return

    url     = sys.argv[1]
    #mixFile = prepareDileptonCollection(url)    
    mixFile='dilepton_summary.pck'
    outURL  = sys.argv[2]
    createMixedFriendTrees(url,mixFile,outURL)



if __name__ == "__main__":
    main()
