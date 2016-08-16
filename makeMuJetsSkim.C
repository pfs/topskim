#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include "LepJetsSkimTree.h"
#include <string>
#include <vector>

#include "ForestTreeHeaders/ForestMuons.h"

const bool isDebug = false;

// Jet and lepton selection
const float jetPtCut  = 25.;
const float jetEtaCut = 3.;
const int   minNJets  = 3;   //Note: you will have to change this to make control distributions for events with less than 4 jets

const float muEtaCut = 2.1;
const float muPtCut  = 15.;

const float muChi2NDFCut   = 10;
const float muInnerD0Cut   = 0.2;
const float muInnerDzCut   = 20.;//0.5;
const int   muMuonHitsCut  = 0;
const int   muStationsCut  = 1;
const int   muTrkLayersCut = 5;
const int   muPixelHitsCut = 0;

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi);

void makeMuJetsSkim(const std::string outFileName = "", const std::string inFileName = "", bool isMC = false)
{
  if(!strcmp(inFileName.c_str(), "")){
    std::cout << "No inputs specified. return" << std::endl;
    return;
  }

  if(isDebug) std::cout << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  skimTree_p = new TTree("skimTree", "skimTree");
  BookTree();

  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
  inFileNames_p->push_back(inFileName);

  TChain *lepTree_p = new TChain("ggHiNtuplizer/EventTree");
  TChain *jetTree_p = new TChain("akCs2PFJetAnalyzer/t");
  TChain *genTree_p = new TChain("HiGenParticleAna/hi");
  TChain *hiTree_p = new TChain("hiEvtAnalyzer/HiTree");
  TChain *hltTree_p = new TChain("hltanalysis/HltTree");
  TChain *pfTree_p = new TChain("pfcandAnalyzerCS/pfTree");
  TChain *skimAnaTree_p = new TChain("skimanalysis/HltTree");
  
  const int nFiles = (int)inFileNames_p->size();

  for(int fileIter = 0; fileIter < nFiles; fileIter++){
    std::cout << "On file: " << fileIter << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
    lepTree_p->Add(inFileNames_p->at(fileIter).c_str());
    jetTree_p->Add(inFileNames_p->at(fileIter).c_str());
    genTree_p->Add(inFileNames_p->at(fileIter).c_str());
    hiTree_p->Add(inFileNames_p->at(fileIter).c_str());
    hltTree_p->Add(inFileNames_p->at(fileIter).c_str());
    pfTree_p->Add(inFileNames_p->at(fileIter).c_str());
    skimAnaTree_p->Add(inFileNames_p->at(fileIter).c_str());
  }
  
  ForestMuons fForestMu;
      
  const int maxJets = 5000;
  int           nref;
  float         jtpt[maxJets];   //[nref]
  float         jteta[maxJets];   //[nref]
  float         jtphi[maxJets];   //[nref]
  float         jtm[maxJets];   //[nref]
  float         discr_csvV1[maxJets]; //[nref]
  float         discr_csvV2[maxJets];
  float         discr_tcHighEff[maxJets];
  float         discr_tcHighPur[maxJets];
  int         refparton_flavorForB[maxJets];

  std::vector<int>     *genpdg = 0;
  std::vector<float>   *genpt = 0;
  std::vector<float>   *geneta = 0;
  std::vector<float>   *genphi = 0;
  std::vector<int>     *genchg = 0;

  //pf particles pfId, pfPt, pfEta, pfPhi
  std::vector<int>           *pfId = 0;
  std::vector<float>         *pfPt = 0;
  std::vector<float>         *pfEta = 0;
  std::vector<float>         *pfPhi = 0;
    
  int trig = 1;

  //event selections
  int phfCoincFilter = 1;
  int HBHENoiseFilterResult = 1;
  int pprimaryVertexFilter = 1;
  int pcollisionEventSelection = 1;

  //lepTree_p->SetBranchStatus("*", 0);
  lepTree_p->SetBranchStatus("mu*", 1);
      
  lepTree_p->SetBranchAddress("muPt", &fForestMu.muPt);
  lepTree_p->SetBranchAddress("muPhi", &fForestMu.muPhi);
  lepTree_p->SetBranchAddress("muEta", &fForestMu.muEta);
  lepTree_p->SetBranchAddress("muCharge", &fForestMu.muCharge);
  lepTree_p->SetBranchAddress("muChi2NDF", &fForestMu.muChi2NDF);
  lepTree_p->SetBranchAddress("muInnerD0", &fForestMu.muInnerD0);
  lepTree_p->SetBranchAddress("muInnerDz", &fForestMu.muInnerDz);
  lepTree_p->SetBranchAddress("muMuonHits", &fForestMu.muMuonHits);
  lepTree_p->SetBranchAddress("muStations", &fForestMu.muStations);
  lepTree_p->SetBranchAddress("muTrkLayers", &fForestMu.muTrkLayers);
  lepTree_p->SetBranchAddress("muPixelHits", &fForestMu.muPixelHits);    
 
  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtm", 1);
  jetTree_p->SetBranchStatus("discr_csvV1", 1);
  jetTree_p->SetBranchStatus("discr_csvV2", 1);
  jetTree_p->SetBranchStatus("discr_tcHighEff", 1);
  jetTree_p->SetBranchStatus("discr_tcHighPur", 1);
  jetTree_p->SetBranchStatus("refparton_flavorForB", 1);
        
  jetTree_p->SetBranchAddress("nref", &nref);
  jetTree_p->SetBranchAddress("jtpt", jtpt);
  jetTree_p->SetBranchAddress("jtphi", jtphi);
  jetTree_p->SetBranchAddress("jteta", jteta);
  jetTree_p->SetBranchAddress("jtm", jtm);
  jetTree_p->SetBranchAddress("discr_csvV1", discr_csvV1);
  jetTree_p->SetBranchAddress("discr_csvV2", discr_csvV2);
  jetTree_p->SetBranchAddress("discr_tcHighEff", discr_tcHighEff);
  jetTree_p->SetBranchAddress("discr_tcHighPur", discr_tcHighPur);
  jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);

  genTree_p->SetBranchStatus("*", 0);
  genTree_p->SetBranchStatus("pdg", 1);
  genTree_p->SetBranchStatus("pt", 1);
  genTree_p->SetBranchStatus("eta", 1);
  genTree_p->SetBranchStatus("phi", 1);
  genTree_p->SetBranchStatus("chg", 1);

  genTree_p->SetBranchAddress("pdg", &genpdg);
  genTree_p->SetBranchAddress("pt", &genpt);
  genTree_p->SetBranchAddress("eta", &geneta);
  genTree_p->SetBranchAddress("phi", &genphi);
  genTree_p->SetBranchAddress("chg", &genchg);
    
  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("vz", 1);
  hiTree_p->SetBranchStatus("weight", 1);
    
  hiTree_p->SetBranchAddress("run", &run_);
  hiTree_p->SetBranchAddress("evt", &evt_);
  hiTree_p->SetBranchAddress("lumi", &lumi_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vz", &vz_);
  hiTree_p->SetBranchAddress("weight", &weight_);

  pfTree_p->SetBranchAddress("pfId", &pfId);
  pfTree_p->SetBranchAddress("pfPt", &pfPt);
  pfTree_p->SetBranchAddress("pfEta", &pfEta);
  pfTree_p->SetBranchAddress("pfPhi", &pfPhi);
    
  hltTree_p->SetBranchStatus("HLT_HIL2Mu15_v2",1);
  hltTree_p->SetBranchAddress("HLT_HIL2Mu15_v2",&trig);

  skimAnaTree_p->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  skimAnaTree_p->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  skimAnaTree_p->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  skimAnaTree_p->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    
  if(isDebug) std::cout << __LINE__ << std::endl;
    
  if(isDebug) std::cout << __LINE__ << std::endl;
    
  int nEntries = (int)lepTree_p->GetEntries();
  //nEntries = 5000;
  //nEntries = 300;
  int entryDiv = ((int)(nEntries/20));
    
  if(isDebug) std::cout << __LINE__ << std::endl;
  

/*Editted by Luuk*/
  //Used to count the percentage of decrease after eacht cut. Are printed at end of the script
  int t1(0), t2(0), t3(0), t4(0), t5(0), t6(0), t7(0), t8(0);
  int m1(0), m2(0), m3(0), m4(0), m5(0), m6(0), m7(0), m8(0), m9(0), m10(0), m11(0), m12(0), m13(0), m14(0);
/**/

  for(int entry = 0; entry < nEntries; entry++){
    if(isDebug) std::cout << __LINE__ << std::endl;
      
    if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    if(isDebug) std::cout << __LINE__ << std::endl;

    hiTree_p->GetEntry(entry);
    lepTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);
    genTree_p->GetEntry(entry);
    hltTree_p->GetEntry(entry);
    pfTree_p->GetEntry(entry);
    skimAnaTree_p->GetEntry(entry);

    if(!isMC && !trig) continue;
    t1++;
    if(!phfCoincFilter) continue;
    t2++;
    if(!HBHENoiseFilterResult) continue;
    t3++;
    if(!pcollisionEventSelection) continue;
    t4++;
    if(!pprimaryVertexFilter) continue;
    t5++;
      
    if(TMath::Abs(vz_) > 15) continue;
    t6++;
      
    if(isDebug) std::cout << __LINE__ << std::endl;

    float tempMuPt_[nLep];
    float tempMuPhi_[nLep];
    float tempMuEta_[nLep];
    int tempMuChg_[nLep];
    float tempMuIso_[nLep];     
    float tempMuInnerDz_[nLep];

    for(int lepIter = 0; lepIter < nLep; lepIter++){
      lepPt_[lepIter] = -999;
      lepPhi_[lepIter] = -999;
      lepEta_[lepIter] = -999;
      lepChg_[lepIter] = -999;
      lepID_[lepIter] = -999;
      lepIso_[lepIter] = -999;
      lepInnerDz_[lepIter] = -999;
    }
       
    for(int lepIter = 0; lepIter < 2; lepIter++){
      tempMuPt_[lepIter] = -999;
      tempMuPhi_[lepIter] = -999;
      tempMuEta_[lepIter] = -999;
      tempMuChg_[lepIter] = -999;
      tempMuIso_[lepIter] = -999;
      tempMuInnerDz_[lepIter] = -999;
    }
       
    for(int ij = 0; ij<nMaxJets; ++ij) {
      jtPt_[ij] = -999.;
      jtEta_[ij] = -999.;
      jtPhi_[ij] = -999.;
      jtM_[ij] = -999.;
      discr_csvV1_[ij] = -999.;
      discr_csvV2_[ij] = -999.;
      discr_tcHighEff_[ij] = -999.;
      discr_tcHighPur_[ij] = -999.;
      refparton_flavorForB_[ij] = -999;
    }

    for(int ijgen = 0; ijgen < nMaxGen; ++ijgen){
      genPdg_[ijgen] = -999;
      genPt_[ijgen] = -999;
      genEta_[ijgen] = -999;
      genPhi_[ijgen] = -999;
      genChg_[ijgen] = -999;
    }
      
    if(isDebug) std::cout << __LINE__ << std::endl;
      
    //Find two leading muons
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {
      if(TMath::Abs(fForestMu.muEta->at(muIter)) > muEtaCut) continue;
      m1++;
      if(fForestMu.muPt->at(muIter) < muPtCut) continue;
      m2++;
  
      if(fForestMu.muChi2NDF->at(muIter) >= muChi2NDFCut) continue;
      m3++;
      if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) > muInnerD0Cut) continue;
      m4++;
      if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) > muInnerDzCut) continue;
      m5++;
      if(fForestMu.muMuonHits->at(muIter) <= muMuonHitsCut) continue;
      m6++;
      if(fForestMu.muStations->at(muIter) <= muStationsCut) continue;
      m7++;
      if(fForestMu.muTrkLayers->at(muIter) <= muTrkLayersCut) continue;
      m8++;
      if(fForestMu.muPixelHits->at(muIter) <= muPixelHitsCut) continue;
      m9++;

      if(fForestMu.muPt->at(muIter) > lepPt_[0]){
        tempMuPt_[1]  = tempMuPt_[0];
        tempMuPhi_[1] = tempMuPhi_[0];
        tempMuEta_[1] = tempMuEta_[0];
        tempMuChg_[1] = tempMuChg_[0];
        tempMuIso_[1] = tempMuIso_[0]; 
        tempMuInnerDz_[1] = tempMuInnerDz_[0];
 
        tempMuPt_[0]  = fForestMu.muPt->at(muIter);
        tempMuPhi_[0] = fForestMu.muPhi->at(muIter);
        tempMuEta_[0] = fForestMu.muEta->at(muIter);
        tempMuChg_[0] = fForestMu.muCharge->at(muIter);
        tempMuIso_[0] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi); //iso;
        tempMuInnerDz_[0] = fForestMu.muInnerDz->at(muIter);
      }
      else if(fForestMu.muPt->at(muIter) > tempMuPt_[1]){
        tempMuPt_[1]  = fForestMu.muPt->at(muIter);
        tempMuPhi_[1] = fForestMu.muPhi->at(muIter);
        tempMuEta_[1] = fForestMu.muEta->at(muIter);
        tempMuChg_[1] = fForestMu.muCharge->at(muIter);
        tempMuIso_[1] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi);//iso;
        tempMuInnerDz_[1] = fForestMu.muInnerDz->at(muIter);
      }
    }

    //store muons in out tree
    int lepIter = 0;
    for(int muIter = 0; muIter < 2; muIter++){

      if(tempMuPt_[muIter]<0.) continue;
      m11++;

      //if(hiBin_<20 && tempMuIso_[muIter]>0.58) continue;
      //else if(hiBin_>=20 && hiBin_<60 && tempMuIso_[muIter]>0.45) continue;
      //else if(hiBin_>=60 && hiBin_<100 && tempMuIso_[muIter]>0.3) continue;
      //else if(hiBin_>=100 && hiBin_<140 && tempMuIso_[muIter]>0.24) continue;
      //else if(hiBin_>=140 && tempMuIso_[muIter]>0.18) continue;
      m12++;

      lepPt_[lepIter] = tempMuPt_[muIter];
      lepPhi_[lepIter] = tempMuPhi_[muIter];
      lepEta_[lepIter] = tempMuEta_[muIter];
      lepChg_[lepIter] = tempMuChg_[muIter];
      lepID_[lepIter] = muID;
      lepIso_[lepIter] = tempMuIso_[muIter];
      lepInnerDz_[lepIter] = tempMuInnerDz_[muIter];
      ++lepIter;
    }
    if(lepIter<1) continue;
    t7++;
    m13 += lepIter;
    nLep_ = lepIter;

    int njets = 0;
    for(int jetIter = 0; jetIter < nref; jetIter++){
      if(jtpt[jetIter]<jetPtCut) continue;
      if(fabs(jteta[jetIter])>jetEtaCut) continue;
      jtPt_[njets]  = jtpt[jetIter];
      jtEta_[njets] = jteta[jetIter];
      jtPhi_[njets] = jtphi[jetIter];
      jtM_[njets]   = jtm[jetIter]; 
      discr_csvV1_[njets] = discr_csvV1[jetIter];
      discr_csvV2_[njets] = discr_csvV2[jetIter];
      discr_tcHighEff_[njets] = discr_tcHighEff[jetIter];
      discr_tcHighPur_[njets] = discr_tcHighPur[jetIter];
      refparton_flavorForB_[njets] = refparton_flavorForB[jetIter];
      ++njets;
    }
    nJt_ = njets;
    if(nJt_<minNJets) continue; //need at least 2 b jets (t->Wb) and 2 light jets (W->qqbar)
    t8++;
    m14 += lepIter;

    std::copy(genpdg->begin(), genpdg->end(), genPdg_);
    std::copy(genpt->begin(), genpt->end(), genPt_);
    std::copy(geneta->begin(), geneta->end(), genEta_);
    std::copy(genphi->begin(), genphi->end(), genPhi_);
    std::copy(genchg->begin(), genchg->end(), genChg_);
    nGen_ = (int)genpdg->size();

    skimTree_p->Fill();
    
  }//entries loop

  
  outFile_p->cd();
  TNamed pathStr1("pathStr1", outFileName.c_str());
  pathStr1.Write("", TObject::kOverwrite);
  TNamed pathStr2("pathStr2", inFileName.c_str());
  pathStr2.Write("", TObject::kOverwrite);
  skimTree_p->Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  cout << "nEntries = " << nEntries << endl;
  cout << "Events: t1 = " << t1 << " , t2 = " << t2 << " , t3 = " << t3 << " , t4 = " << t4 << " , t5 = " << t5 << " , t6 = " << t6 << " , t7 = " << t7 << " , t8 = " << t8 << endl;
  cout << "Muons: m1 = " << m1 << " , m2 = " << m2 << " , m3 = " << m3 << " , m4 = " << m4 << " , m5 = " << m5 << " , m6 = " << m6 << " , m7 = " << m7 << " , m8 = " << m8 << " , m9 = " << m9 << " , m10 = X" << " , m11 = " << m11 << " , m12 = " << m12 << " , m13 = " << m13 << " , m14 = " << m14 << endl;

  return;
}

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi) {
  //calculate lepton isolation from pf candidates.
  //Isolation cone R=0.3
  //excluding pf candidates at a distance less than 0.03 from lepton
  
  double conePt = 0.;
  for(unsigned int i = 0; i<pfPt->size(); ++i) {

    double deltaR = sqrt(pow(acos(cos(lepPhi-pfPhi->at(i))),2)+pow(lepEta-pfEta->at(i),2));

    if(deltaR<0.03 || deltaR>0.3) continue;

    conePt+=pfPt->at(i);
  }
  double relIso = conePt;
  if(lepPt>0.) relIso = conePt/lepPt;
  
  return relIso;
}
