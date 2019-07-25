#ifndef ForestLeptons_h
#define ForestLeptons_h

#include <iostream>
#include <vector>

class ForestLeptons {
public :
 ForestLeptons(TChain *t) 
   : 
    eleChargeConsistent(0),
    eleEn(0),
    eleDzErr(0),
    eleTrkPt(0),
    eleTrkEta(0),
    eleTrkPhi(0),
    eleTrkCharge(0),
    eleTrkChi2(0),
    eleTrkNdof(0),
    eleTrkNormalizedChi2(0),
    eleTrkValidHits(0),
    eleTrkLayers(0),
    elePt(0),
    eleEta(0),
    elePhi(0),
    eleSCEn(0),
    eleESEn(0),
    eleSCEta(0),
    eleSCPhi(0),
    eleSCRawEn(0),
    eleSCEtaWidth(0),
    eleSCPhiWidth(0),
    eleHoverE(0),
    eleHoverEBc(0),
    eleEoverP(0),
    eleEoverPInv(0),
    eleBrem(0),
    eledEtaAtVtx(0),
    eledEtaSeedAtVtx(0),
    eledPhiAtVtx(0),
    eleSigmaIEtaIEta(0),
    eleSigmaIEtaIEta_2012(0),
    eleSigmaIPhiIPhi(0),
    eleMissHits(0),
    eleESEffSigmaRR(0),
    elePFChIso(0),
    elePFPhoIso(0),
    elePFNeuIso(0),
    elePFPUIso(0),
    elePFChIso03(0),
    elePFPhoIso03(0),
    elePFNeuIso03(0),
    elePFChIso04(0),
    elePFPhoIso04(0),
    elePFNeuIso04(0),
    eleD0(0),
    eleIP3D(0),
    eleCharge(0),
    eleD0Err(0),
    eleBC2Eta(0),
    muPt(0),
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

       t->SetBranchStatus("ele*", 1);
       t->SetBranchAddress("elePt", &elePt);
       t->SetBranchAddress("elePhi", &elePhi);
       t->SetBranchAddress("eleEta", &eleEta);
       t->SetBranchAddress("eleCharge", &eleCharge);
       t->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta);
       t->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx);
       t->SetBranchAddress("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
       t->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx);
       t->SetBranchAddress("eleHoverE", &eleHoverE);
       t->SetBranchAddress("eleHoverEBc", &eleHoverEBc);
       t->SetBranchAddress("eleD0", &eleD0);
       t->SetBranchAddress("eleIP3D", &eleIP3D);
       t->SetBranchAddress("eleD0Err", &eleD0Err);
       t->SetBranchAddress("eleDz", &eleDz);
       t->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
       t->SetBranchAddress("eleSCPhi", &eleSCPhi);
       t->SetBranchAddress("eleSCEta", &eleSCEta);
       t->SetBranchAddress("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
       t->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi);
       t->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth);
       t->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth);
       t->SetBranchAddress("eleMissHits", &eleMissHits);
       t->SetBranchAddress("elePFChIso", &elePFChIso);
       t->SetBranchAddress("elePFPhoIso", &elePFPhoIso);
       t->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
       t->SetBranchAddress("elePFChIso03", &elePFChIso03);
       t->SetBranchAddress("elePFPhoIso03", &elePFPhoIso03);
       t->SetBranchAddress("elePFNeuIso03", &elePFNeuIso03);
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
   ~ForestLeptons(){};

   // RecoElectron info
   Int_t           nEle;
   std::vector<Int_t>     *eleChargeConsistent;
   std::vector<Float_t>   *eleEn;
   std::vector<Float_t>   *eleDzErr;
   std::vector<Float_t>   *eleTrkPt;
   std::vector<Float_t>   *eleTrkEta;
   std::vector<Float_t>   *eleTrkPhi;
   std::vector<Int_t>     *eleTrkCharge;
   std::vector<Float_t>   *eleTrkChi2;
   std::vector<Float_t>   *eleTrkNdof;
   std::vector<Float_t>   *eleTrkNormalizedChi2;
   std::vector<Int_t>     *eleTrkValidHits;
   std::vector<Int_t>     *eleTrkLayers;
   std::vector<Float_t>   *elePt;
   std::vector<Float_t>   *eleEta;
   std::vector<Float_t>   *elePhi;
   std::vector<Float_t>   *eleSCEn;
   std::vector<Float_t>   *eleESEn;
   std::vector<Float_t>   *eleSCEta;
   std::vector<Float_t>   *eleSCPhi;
   std::vector<Float_t>   *eleSCRawEn;
   std::vector<Float_t>   *eleSCEtaWidth;
   std::vector<Float_t>   *eleSCPhiWidth;
   std::vector<Float_t>   *eleHoverE;
   std::vector<Float_t>   *eleHoverEBc;
   std::vector<Float_t>   *eleEoverP;
   std::vector<Float_t>   *eleEoverPInv;
   std::vector<Float_t>   *eleBrem;
   std::vector<Float_t>   *eledEtaAtVtx;
   std::vector<Float_t>   *eledEtaSeedAtVtx;
   std::vector<Float_t>   *eledPhiAtVtx;
   std::vector<Float_t>   *eleSigmaIEtaIEta;
   std::vector<Float_t>   *eleSigmaIEtaIEta_2012;
   std::vector<Float_t>   *eleSigmaIPhiIPhi;
   std::vector<Int_t>     *eleMissHits;
   std::vector<Float_t>   *eleESEffSigmaRR;
   std::vector<Float_t>   *elePFChIso;
   std::vector<Float_t>   *elePFPhoIso;
   std::vector<Float_t>   *elePFNeuIso;
   std::vector<Float_t>   *elePFPUIso;
   std::vector<Float_t>   *elePFChIso03;
   std::vector<Float_t>   *elePFPhoIso03;
   std::vector<Float_t>   *elePFNeuIso03;
   std::vector<Float_t>   *elePFChIso04;
   std::vector<Float_t>   *elePFPhoIso04;
   std::vector<Float_t>   *elePFNeuIso04;
   std::vector<Float_t>   *eleD0;
   std::vector<Float_t>   *eleIP3D;
   std::vector<Int_t>     *eleCharge;
   std::vector<Float_t>   *eleD0Err;
   std::vector<Float_t>   *eleDz;
   std::vector<Float_t>   *eleBC2Eta;

   // RecoMuon info
   std::vector<Float_t>   *muPt;
   std::vector<Float_t>   *muEta;
   std::vector<Float_t>   *muPhi;
   std::vector<Int_t>     *muCharge;
   std::vector<Int_t>     *muType;
   std::vector<Int_t>     *muIsGood;
   std::vector<Float_t>   *muD0;
   std::vector<Float_t>   *muDz;
   std::vector<Float_t>   *muChi2NDF;
   std::vector<Float_t>   *muInnerD0;
   std::vector<Float_t>   *muInnerDz;
   std::vector<Int_t>     *muTrkLayers;
   std::vector<Int_t>     *muPixelLayers;
   std::vector<Int_t>     *muPixelHits;
   std::vector<Int_t>     *muMuonHits;
   std::vector<Int_t>     *muTrkQuality;
   std::vector<Int_t>     *muStations;
   std::vector<Float_t>   *muIsoTrk;
   std::vector<Float_t>   *muPFChIso;
   std::vector<Float_t>   *muPFPhoIso;
   std::vector<Float_t>   *muPFNeuIso;
   std::vector<Float_t>   *muPFPUIso;
};
#endif
