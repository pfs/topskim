#ifndef ForestMuons_h
#define ForestMuons_h

#define maxForestMuons 10

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestMuons {
public :
 ForestMuons(TChain *t) : muPt(0),
    muEta(0),
    muPhi(0),
    muCharge(0),
    muType(0),
    muIsGood(0),
    muD0(0),
    muDz(0),
    muChi2NDF(0),
    muInnerD0(0),
    muInnerDz(0),
    muTrkLayers(0),
    muPixelLayers(0),
    muPixelHits(0),
    muMuonHits(0),
    muTrkQuality(0),
    muStations(0),
    muIsoTrk(0),
    muPFChIso(0),
    muPFPhoIso(0),
    muPFNeuIso(0),
    muPFPUIso(0)
      {
        t->SetBranchStatus("mu*", 1);
        t->SetBranchAddress("muPt", &muPt);
        t->SetBranchAddress("muPhi", &muPhi);
        t->SetBranchAddress("muEta", &muEta);
        t->SetBranchAddress("muCharge", &muCharge);
        t->SetBranchAddress("muType", &muType);
        t->SetBranchAddress("muIsGood", &muIsGood);
        t->SetBranchAddress("muD0", &muD0);
        t->SetBranchAddress("muDz", &muDz);
        t->SetBranchAddress("muChi2NDF", &muChi2NDF);
        t->SetBranchAddress("muInnerD0", &muInnerD0);
        t->SetBranchAddress("muInnerDz", &muInnerDz);
        t->SetBranchAddress("muMuonHits", &muMuonHits);
        t->SetBranchAddress("muStations", &muStations);
        t->SetBranchAddress("muTrkLayers", &muTrkLayers);
        t->SetBranchAddress("muPixelHits", &muPixelHits);    
        t->SetBranchAddress("muMuonHits", &muMuonHits);
        t->SetBranchAddress("muTrkQuality", &muTrkQuality);
        t->SetBranchAddress("muStations", &muStations);
        t->SetBranchAddress("muIsoTrk", &muIsoTrk);
        t->SetBranchAddress("muPFChIso", &muPFChIso);
        t->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
        t->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
        t->SetBranchAddress("muPFPUIso", &muPFPUIso);
      };

  ~ForestMuons(){};

   // RecoMuon info
   std::vector<float>   *muPt;
   std::vector<float>   *muEta;
   std::vector<float>   *muPhi;
   std::vector<int>     *muCharge;
   std::vector<int>     *muType;
   std::vector<int>     *muIsGood;
   std::vector<float>   *muD0;
   std::vector<float>   *muDz;
   std::vector<float>   *muChi2NDF;
   std::vector<float>   *muInnerD0;
   std::vector<float>   *muInnerDz;
   std::vector<int>     *muTrkLayers;
   std::vector<int>     *muPixelLayers;
   std::vector<int>     *muPixelHits;
   std::vector<int>     *muMuonHits;
   std::vector<int>     *muTrkQuality;
   std::vector<int>     *muStations;
   std::vector<float>   *muIsoTrk;
   std::vector<float>   *muPFChIso;
   std::vector<float>   *muPFPhoIso;
   std::vector<float>   *muPFNeuIso;
   std::vector<float>   *muPFPUIso;
 };
#endif
