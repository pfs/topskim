#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TSystem.h"

#include <string>
#include <vector>

#include "HeavyIonsAnalysis/topskim/include/ForestHiTree.h"
#include "HeavyIonsAnalysis/topskim/include/ForestElectrons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestMuons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestPFCands.h"
#include "HeavyIonsAnalysis/topskim/include/ForestJets.h"
#include "HeavyIonsAnalysis/topskim/include/LumiRun.h"
#include "HeavyIonsAnalysis/topskim/include/HistTool.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

const bool isDebug = true;

const float jetPtCut  = 30.;
const float jetEtaCut = 2.4;
const float lepPtCut  = 20.;
const float lepEtaCut = 2.1;
//see https://indico.cern.ch/event/803679/contributions/3342407/attachments/1808912/2953435/egm-minipog-190308.pdf
const float eeScaleShift = 6.8182E-2/5.9097E-2;
const int firstEEScaleShiftRun = 327402; 
const float barrelEndcapEta[2]={1.4442,1.5660};
const float csvWP = 0.8838;

using namespace std;
using namespace fastjet;


//note that it is constexpr function
template <typename T>
constexpr typename std::underlying_type<T>::type integral(T value) 
{
  return static_cast<typename std::underlying_type<T>::type>(value);
}

// index, ntks in svtx, m svtx, csv
enum class JetInfo{
                    index=0, ntks=1, msvtx=2, csvV2=3,
                    phi=4, pu=5, mass=6,
                    trackMax=7,   trackSum=8,   trackN=9, trackHardSum=10, trackHardN=11, 
	            chargedMax=12, chargedSum=13, chargedN=14, chargedHardSum=15, chargedHardN=16, 
	            photonMax=17, photonSum=18, photonHardSum=19, photonHardN=20,
	            neutralMax=21,neutralSum=22,neutralN=23,
		    hcalSum=24,   ecalSum=25, 
		    eMax=26,      eSum=27, eN=28, muMax=29,muSum=30, muN=31,
		    PfCHF=32, PfNHF=33, PfCEF=34, PfNEF=35, PfMUF=36, PfCHM=37, PfNHM=38, PfCEM=39, PfNEM=40
		    };
		    
typedef std::tuple<int, int, float, float,
		   float, float, float,
		   float, float, int, float, int,
		   float, float, int, float, int,
		   float, float, float, int,
		   float, float, int, 
		   float, float, 
		   float, float, int, float, float, int,
		   float, float, float, float, float, int, int, int, int
		   > JetInfo_t;

static bool orderByBtagInfo(const JetInfo_t &a, const JetInfo_t &b)
{
  //int ntks_a(std::get<1>(a)), ntks_b(std::get<1>(b));
  //if(ntks_a>ntks_b) return true;

  float csv_a(std::get<integral(JetInfo::csvV2)>(a)), csv_b(std::get<integral(JetInfo::csvV2)>(b));
  if(csv_a>csv_b) return true;
  return false;
}


//
int main(int argc, char* argv[])
{
  bool blind(true);
  TString inURL,outURL;
  bool isMC(false),isPP(false);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--in")!=string::npos && i+1<argc)       { inURL=TString(argv[i+1]); i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc) { outURL=TString(argv[i+1]); i++;}
    else if(arg.find("--mc")!=string::npos)              { isMC=true;  }
    else if(arg.find("--pp")!=string::npos)              { isPP=true;  }
  }

  bool isSingleMuPD( !isMC && inURL.Contains("SkimMuons"));
  bool isSingleElePD( !isMC && inURL.Contains("SkimElectrons"));
  LumiRun lumiTool;

  if(isPP)
    cout << "Treating as a pp collision file" << endl;
  if(isMC)
    cout << "Treating as a MC file" << endl;

  //book some histograms
  HistTool ht;

  if(!isMC) ht.addHist("ratevsrun",lumiTool.getLumiMonitor());

  //generic histograms
  for(int i=0; i<2; i++) {
    TString pf(Form("l%d",i+1));
    ht.addHist(pf+"pt",        new TH1F(pf+"pt",       ";Lepton transverse momentum [GeV];Events",20,20,200));
    ht.addHist(pf+"eta",       new TH1F(pf+"eta",      ";Lepton pseudo-rapidity;Events",40,0,5.0));
    ht.addHist(pf+"chreliso",  new TH1F(pf+"chreliso", ";Relative PF charged isolation;Leptons",50,0,5.0));
    ht.addHist(pf+"phoreliso", new TH1F(pf+"phoreliso",";Relative PF photon isolation;Leptons",50,0,5.0));
    ht.addHist(pf+"neureliso", new TH1F(pf+"neureliso",";Relative PF neutral hadron isolation;Leptons",50,0,5.0));
    ht.addHist(pf+"chrelisovscen",  new TH2F(pf+"chrelisovscen", ";Relative PF charged isolation;Centrality bin;Leptons",20,0,2.0,5,0,100));
    ht.addHist(pf+"phorelisovscen", new TH2F(pf+"phorelisovscen",";Relative PF photon isolation;Centrality bin;Leptons",20,0,1.0,5,0,100));
    ht.addHist(pf+"neurelisovscen", new TH2F(pf+"neurelisovscen",";Relative PF neutral hadron isolation;Centrality bin;Leptons",20,0,1.0,5,0,100));
  }
  ht.addHist("mll",      new TH1F("mll",      ";Dilepton invariant mass [GeV];Events",20,20,200));
  ht.addHist("ptll",     new TH1F("ptll",     ";Dilepton transverse momentum [GeV];Events",20,0,200));
  ht.addHist("dphill",   new TH1F("dphill",   ";#Delta#phi(l,l');Events",20,0,3.15));
  ht.addHist("detall",   new TH1F("detall",   ";#Delta#eta(l,l');Events",20,0,5));
  ht.addHist("chrho",    new TH1F("chrho",    ";#rho_{ch};Events",25,0,25));
  for(size_t i=0; i<2; i++) {
    TString pf(i==0 ? "tk" : "pf");
    ht.addHist("n"+pf+"jets",    new TH1F("n"+pf+"jets",    ";Jet multiplicity;Events",25,0,25));
    ht.addHist("n"+pf+"bjets",   new TH1F("n"+pf+"bjets",   ";b-jet multiplicity;Events",25,0,25));    
    ht.addHist("n"+pf+"svtx",    new TH1F("n"+pf+"svtx",    ";Secondary vertex multiplicity;Events",15,0,15));
    for(size_t j=1; j<=2; j++){
      TString ppf(j==1 ? "1" : "2");
      ht.addHist(pf+ppf+"jbalance",    new TH1F(pf+ppf+"jbalance", ";R = p_{T}(ll)/p_{T}(j);Events",150,0,5.5));
      ht.addHist(pf+ppf+"jpt",         new TH1F(pf+ppf+"jpt",      ";Jet transverse momentum [GeV];Events",20,30,200));
      ht.addHist(pf+ppf+"jeta",        new TH1F(pf+ppf+"jeta",     ";Jet pseudo-rapidity;Events",40,0,5.0));
      ht.addHist(pf+ppf+"jphi",        new TH1F(pf+ppf+"jphi",     ";Jet azimuthal angle;Events",40,-3.15,3.15));
      ht.addHist(pf+ppf+"jpu",         new TH1F(pf+ppf+"jpu",      ";Jet pu [GeV];Events",20,0,100));
      ht.addHist(pf+ppf+"jm",          new TH1F(pf+ppf+"jm",       ";Jet mass [GeV];Events",20,0,20));
      ht.addHist(pf+ppf+"jPfCHF",      new TH1F(pf+ppf+"jPfCHF",    ";Jet PfCHF;Events",25,0,1.));
      ht.addHist(pf+ppf+"jPfNHF",      new TH1F(pf+ppf+"jPfNHF",    ";Jet PfNHF;Events",25,0,1.));
      ht.addHist(pf+ppf+"jPfCEF",      new TH1F(pf+ppf+"jPfCEF",    ";Jet PfCEF;Events",25,0,1.));
      ht.addHist(pf+ppf+"jPfNEF",      new TH1F(pf+ppf+"jPfNEF",    ";Jet PfNEF;Events",25,0,1.));
      ht.addHist(pf+ppf+"jPfMUF",      new TH1F(pf+ppf+"jPfMUF",    ";Jet PfMUF;Events",25,0,1.));
      ht.addHist(pf+ppf+"jPfCHM",      new TH1F(pf+ppf+"jPfCHM",    ";Jet PfCHM;Events",25,0,50));
      ht.addHist(pf+ppf+"jPfNHM",      new TH1F(pf+ppf+"jPfNHM",    ";Jet PfNHM;Events",10,0,10));
      ht.addHist(pf+ppf+"jPfCEM",      new TH1F(pf+ppf+"jPfCEM",    ";Jet PfCEM;Events",10,0,10));
      ht.addHist(pf+ppf+"jPfNEM",      new TH1F(pf+ppf+"jPfNEM",    ";Jet PfNEM,;Events",20,0,20));
      ht.addHist(pf+ppf+"jPfMUM",      new TH1F(pf+ppf+"jPfMUM",    ";Jet PfMUM;Events",10,0,10));
      ht.addHist(pf+ppf+"jtrackMax",   new TH1F(pf+ppf+"jtrackMax", ";Jet trackMax [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"jtrackSum",   new TH1F(pf+ppf+"jtrackSum", ";Jet trackSum [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"jtrackN",     new TH1F(pf+ppf+"jtrackN",   ";Jet trackN;Events",25,0,200));
      ht.addHist(pf+ppf+"jtrackHardSum",     new TH1F(pf+ppf+"jtrackHardSum",    ";Jet trackHardSum [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"jtrackHardN",       new TH1F(pf+ppf+"jtrackHardN",    ";Jet trackHardN;Events",10,0,10.));
      ht.addHist(pf+ppf+"chargedMax",        new TH1F(pf+ppf+"jchargedMax",    ";Jet chargedMax [GeV];Events",25,0,200.));
      ht.addHist(pf+ppf+"chargedSum",        new TH1F(pf+ppf+"jchargedSum",    ";Jet chargedSum [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"chargedN",          new TH1F(pf+ppf+"jchargedN",    ";Jet chargedN;Events",25,0,200));
      ht.addHist(pf+ppf+"chargedHardSum",    new TH1F(pf+ppf+"jchargedHardSum",    ";Jet chargedHardSum [GeV],;Events",20,0,200));
      ht.addHist(pf+ppf+"chargedHardN",      new TH1F(pf+ppf+"jchargedHardN",    ";Jet chargedHardN;Events",10,0,10));
      ht.addHist(pf+ppf+"photonMax",         new TH1F(pf+ppf+"jphotonMax",    ";Jet photonMax [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"photonSum",         new TH1F(pf+ppf+"jphotonSum",    ";Jet photonSum [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"photonN",           new TH1F(pf+ppf+"jphotonN",    ";Jet photonN;Events",25,0,200.));
      ht.addHist(pf+ppf+"photonHardSum",     new TH1F(pf+ppf+"jphotonHardSum",    ";Jet photonHardSum [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"photonHardN",       new TH1F(pf+ppf+"jphotonHardN",    ";Jet photonHardN;Events",10,0,10.));
      ht.addHist(pf+ppf+"neutralMax",        new TH1F(pf+ppf+"jneutralMax",    ";Jet neutralMax [GeV];Events",25,0,200.));
      ht.addHist(pf+ppf+"neutralSum",        new TH1F(pf+ppf+"jneutralSum",    ";Jet neutralSum [GeV];Events",25,0,200));
      ht.addHist(pf+ppf+"neutralN",          new TH1F(pf+ppf+"jneutralN",    ";Jet neutralN;Events",10,0,10));
      ht.addHist(pf+ppf+"hcalSum",           new TH1F(pf+ppf+"jhcalSum",    ";Jet hcalSum [GeV];Events",20,0,200));
      ht.addHist(pf+ppf+"ecalSum",           new TH1F(pf+ppf+"jecalSum",    ";Jet ecalSum;Events",20,0,200));
      ht.addHist(pf+ppf+"eMax",              new TH1F(pf+ppf+"jeMax",    ";Jet eMax [GeV];Events",25,0,200.));
      ht.addHist(pf+ppf+"eSum",              new TH1F(pf+ppf+"jeSum",    ";Jet eSum [GeV];Events",25,0,200.));
      ht.addHist(pf+ppf+"eN",                new TH1F(pf+ppf+"jeN",    ";Jet eN;Events",10,0,10.));
      ht.addHist(pf+ppf+"muMax",             new TH1F(pf+ppf+"jmuMax",    ";Jet muMax [GeV];Events",25,0,200.));
      ht.addHist(pf+ppf+"muSum",             new TH1F(pf+ppf+"jmuSum",    ";Jet muSum [GeV];Events",10,0,10));
      ht.addHist(pf+ppf+"muN",               new TH1F(pf+ppf+"jmuN",    ";Jet muN;Events",10,0,10));
      ht.addHist(pf+ppf+"jsvtxm",            new TH1F(pf+ppf+"jsvtxm",   ";Secondary vertex mass;Events",25,0,6));
      ht.addHist(pf+ppf+"jsvtxntk",          new TH1F(pf+ppf+"jsvtxntk", ";Secondary vertex track multiplicity;Events",5,0,5));
      ht.addHist(pf+ppf+"jcsv",              new TH1F(pf+ppf+"jcsv",     ";CSVv2;Events",25,0,1));
      
    }
  }

  //configure leptons
  TChain *lepTree_p     = new TChain(isPP ? "ggHiNtuplizer/EventTree" : "ggHiNtuplizerGED/EventTree");
  lepTree_p->Add(inURL);
  ForestMuons fForestMu(lepTree_p);  
  ForestElectrons fForestEle(lepTree_p);

  //configure PF cands
  TChain *pfCandTree_p  = new TChain("pfcandAnalyzer/pfTree");
  pfCandTree_p->Add(inURL);
  ForestPFCands fForestPF(pfCandTree_p);

  //configure jets
  TChain *jetTree_p     = new TChain(isPP ? "ak4PFJetAnalyzer/t" : "akPu4PFJetAnalyzer/t");
  jetTree_p->Add(inURL);
  ForestJets fForestJets(jetTree_p);

  //global variables
  TChain *hiTree_p      = new TChain("hiEvtAnalyzer/HiTree");
  hiTree_p->Add(inURL);
  HiTree fForestTree(hiTree_p);

  //trigger
  TChain *hltTree_p     = new TChain("hltanalysis/HltTree");
  hltTree_p->Add(inURL);
  int etrig(0),mtrig(0);
  if(isPP){
    hltTree_p->SetBranchStatus("HLT_HIL3Mu20_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu20_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20_WPLoose_Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20_WPLoose_Gsf_v1",&etrig);
  }else{
    hltTree_p->SetBranchStatus("HLT_HIL3Mu15_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu15_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20Gsf_v1",&etrig);    
  }
    
  Float_t wgtSum(0);
  int nEntries = (int)lepTree_p->GetEntries();  
  int entryDiv = ((int)(nEntries/20));    
  cout << inURL << " has " << nEntries << "events to process" << endl;
  for(int entry = 0; entry < nEntries; entry++){
    
    if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    //cout << "debug 1" << endl;
    lepTree_p->GetEntry(entry);
    //cout << "debug 2" << endl;
    //pfCandTree_p -> Print();
    pfCandTree_p->GetEntry(entry);
    //cout << "debug 3" << endl;
    jetTree_p->GetEntry(entry);    
    //cout << "debug 4" << endl;
    hltTree_p->GetEntry(entry);
    //cout << "debug 5" << endl;
    hiTree_p->GetEntry(entry);

    wgtSum += fForestTree.weight;

    //first of all require a trigger
    int trig=etrig+mtrig;
    if(trig==0) continue;

    //apply global filters
    if(!isPP){
      if(TMath::Abs(fForestTree.vz) > 15) continue;
    }

    float cenBin=0;
    if(!isMC){
      cenBin=0.5*fForestTree.hiBin;
      Int_t runBin=lumiTool.getRunBin(fForestTree.run);
      Float_t lumi=lumiTool.getLumi(fForestTree.run);
      if(lumi>0.){
        if(etrig>0) ht.fill("ratevsrun",runBin,1./lumi,"e");
        if(mtrig>0) ht.fill("ratevsrun",runBin,1./lumi,"m");
      }
    }


    //select muons
    std::vector<int> muIdx,noIdMuIdx;
    std::vector<TLorentzVector> muP4;
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),0.1057);
      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(p4.Pt() < lepPtCut) continue;

      noIdMuIdx.push_back(muIter);

      //id (Tight muon requirements)
      int type=fForestMu.muType->at(muIter);
      bool isGlobal( ((type>>1)&0x1) );
      if(!isGlobal) continue;
      bool isPF( ((type>>5)&0x1) );
      if(!isPF) continue;
      bool isGlobalMuonPromptTight(fForestMu.muChi2NDF->at(muIter)<10. && fForestMu.muMuonHits->at(muIter)>0);
      if(!isGlobalMuonPromptTight) continue;
      if(fForestMu.muStations->at(muIter)<=1) continue;
      if(fForestMu.muTrkLayers->at(muIter) <= 5) continue;
      if(fForestMu.muPixelHits->at(muIter) == 0) continue;
      if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) >=0.2 ) continue;
      if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) >=0.5) continue;

      //selected a good muon
      muIdx.push_back(muIter);
      muP4.push_back(p4);
    }
    
    //select electrons
    //cf. https://twiki.cern.ch/twiki/pub/CMS/HiHighPt2019/HIN_electrons2018_followUp.pdf
    std::vector<int> eleIdx,noIdEleIdx; 
    std::vector<TLorentzVector> eP4;
    bool allEleInEB(true);
    for(unsigned int eleIter = 0; eleIter < fForestEle.elePt->size(); ++eleIter) {
      
      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),0.000511);

      //apply ad-hoc shift for endcap electrons if needed
      if(!isPP && fForestTree.run<=firstEEScaleShiftRun && TMath::Abs(p4.Eta())>barrelEndcapEta[1])
        p4 *=eeScaleShift;         

      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(TMath::Abs(p4.Eta()) > barrelEndcapEta[0] && TMath::Abs(p4.Eta()) < barrelEndcapEta[1] ) continue;
      if(p4.Pt() < lepPtCut) continue;	      

      noIdEleIdx.push_back(eleIter);

      //electron id
      if(fForestEle.eleMissHits->at(eleIter)>1) continue;
      if(fForestEle.eleEoverPInv->at(eleIter)>=0.3) continue;
      if(fForestEle.eleHoverE->at(eleIter)>=0.2) continue;
      if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter))>=0.1) continue;
      if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter))>=0.2) continue;
      if(fForestEle.eleSigmaIEtaIEta->at(eleIter)>=0.05) continue;
      if(TMath::Abs(fForestEle.eleD0->at(eleIter))>=0.1) continue;
      if(TMath::Abs(fForestEle.eleDz->at(eleIter))>=0.5) continue;

      //selected electron
      eleIdx.push_back(eleIter);
      eP4.push_back(p4);
      if(TMath::Abs(p4.Eta()) > barrelEndcapEta[0]) allEleInEB=false;
    }

    int nLep=muIdx.size()+eleIdx.size();
    if(nLep<2) continue;
  
    std::vector<TLorentzVector> selLeptons;
    int dilCode(0);
    int charge(0);    
    TLorentzVector ll;
    vector< std::tuple<float,float,float> > liso;
    if(muP4.size()>1 && mtrig>0) {

      //muon final states from the muon PD only
      if(!isMC && !isSingleMuPD) continue;

      dilCode=13*13;
      ll=muP4[0]+muP4[1];
      selLeptons.push_back(muP4[0]);
      selLeptons.push_back(muP4[1]);
      charge=fForestMu.muCharge->at(muIdx[0])*fForestMu.muCharge->at(muIdx[1]);
      for(size_t i=0; i<2; i++){
        float chIso( fForestMu.muPFChIso->at(muIdx[i]) );
        float phoIso( fForestMu.muPFPhoIso->at(muIdx[i]) );
        float neutIso( fForestMu.muPFNeuIso->at(muIdx[i]) );        
        liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
      }
    }
    else if(muP4.size()>0 && eP4.size()>0 && (etrig>0 || mtrig>0)) {

      //simultaneous triggering only in the single muon PD
      //to avoid double counting
      if(!isMC && etrig>0 && mtrig>0 && isSingleElePD) continue;

      dilCode=11*13;
      ll=muP4[0]+eP4[0];
      selLeptons.push_back(muP4[0]);
      selLeptons.push_back(eP4[0]);
      charge=fForestMu.muCharge->at(muIdx[0])*fForestEle.eleCharge->at(eleIdx[0]);      
      float chIso( fForestMu.muPFChIso->at(muIdx[0]) );
      float phoIso( fForestMu.muPFPhoIso->at(muIdx[0]) );
      float neutIso( fForestMu.muPFNeuIso->at(muIdx[0]) );        
      liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
      chIso=fForestEle.elePFChIso03->at(eleIdx[0]);
      phoIso=fForestEle.elePFPhoIso03->at(eleIdx[0]);
      neutIso=fForestEle.elePFNeuIso03->at(eleIdx[0]);        
      liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
    }
    else if(eP4.size()>1 && etrig>0) {

      //ee final states from the electron PD only 
      if(!isMC && !isSingleElePD) continue;

      dilCode=11*11;
      ll=eP4[0]+eP4[1];
      selLeptons.push_back(eP4[0]);
      selLeptons.push_back(eP4[1]);
      charge=fForestEle.eleCharge->at(eleIdx[0])*fForestEle.eleCharge->at(eleIdx[1]);
      for(size_t i=0; i<2; i++){
        float chIso=fForestEle.elePFChIso03->at(eleIdx[i]);
        float phoIso=fForestEle.elePFPhoIso03->at(eleIdx[i]);
        float neutIso=fForestEle.elePFNeuIso03->at(eleIdx[i]);        
        liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
      }
    }else{
      continue;
    }

    if(ll.M()<20) continue;
    bool isZ( dilCode!=11*13 && fabs(ll.M()-91)<15 );

    if(blind) {
      if(!isMC && !isZ && charge<0 && fForestTree.run>=326887) continue;
    }
    
    TString dilCat("em");
    if(dilCode==11*11) { dilCat=isZ ? "zee" : "ee"; }
    if(dilCode==13*13) { dilCat=isZ ? "zmm" : "mm"; }    
    if(charge>0) dilCat="ss"+dilCat;

    //build track jets from PF candidates
    //cross-clean with respect to the selected leptons
    //require at least 2 constituents
    std::vector<TLorentzVector> tkJetsP4;
    std::vector<PseudoJet> pseudoParticles;
    TLorentzVector p4(0,0,0,0);
    for(size_t ipf=0; ipf<fForestPF.pfId->size(); ipf++) {
      int id(abs(fForestPF.pfId->at(ipf)));

      //pass all neutrals
      if(id==22 || id==130 || id==2112 || id==1 || id==2) continue;      

      //treat all as pions...
      p4.SetPtEtaPhiM(fForestPF.pfPt->at(ipf),fForestPF.pfEta->at(ipf),fForestPF.pfPhi->at(ipf),0.13957);

      //some basic kinematic cuts
      if(p4.Pt()<0.5) continue;
      if(fabs(p4.Eta())<2.5) continue;

      PseudoJet ip=PseudoJet(p4.Px(),p4.Py(),p4.Pz(),p4.E());
      ip.set_user_index( ipf );
      pseudoParticles.push_back( ip );
    }
    JetDefinition jet_def(antikt_algorithm, 0.4);
    ClusterSequence cs(pseudoParticles, jet_def);
    Selector sel_rapmax = SelectorAbsRapMax(2.4);
    JetDefinition jet_def_for_rho(kt_algorithm,0.5);
    AreaDefinition area_def(active_area,GhostedAreaSpec(2.4+1.));
    JetMedianBackgroundEstimator bge(sel_rapmax, jet_def_for_rho, area_def);
    bge.set_particles(pseudoParticles);
    float tkrho=bge.rho();

    std::vector<PseudoJet> tkjets = sorted_by_pt(cs.inclusive_jets());
    for(auto j : tkjets) {
      if(j.constituents().size()<2) continue;
      TLorentzVector p4(j.px(),j.py(),j.pz(),j.e());
      if(p4.DeltaR(selLeptons[0])<0.4 || p4.DeltaR(selLeptons[1])<0.4) continue;
      if(fabs(p4.Eta())<2.4) continue;
      tkJetsP4.push_back(p4);
    }
   
    //b-tag jet the track jets by matching in deltaR to PF jets
    std::vector<JetInfo_t> matchedJetsIdx,pfJetsIdx;
    std::vector<TLorentzVector> pfJetsP4;
    int npfjets(0),npfbjets(0); 
    bool allPFBInEB(true),hasAwayPFJet(true);
    for(int jetIter = 0; jetIter < fForestJets.nref; jetIter++){

      //at least two tracks
      if(fForestJets.trackN[jetIter]<2) continue;

      TLorentzVector jp4(0,0,0,0);
      jp4.SetPtEtaPhiM( fForestJets.jtpt[jetIter],fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);

      float csvVal=fForestJets.discr_csvV2[jetIter];
      int nsvtxTk=fForestJets.svtxntrk[jetIter];
      float msvtx=fForestJets.svtxm[jetIter];
      float phi=fForestJets.jtphi[jetIter];
      float pu=fForestJets.jtpu[jetIter];
      float mass=fForestJets.jtm[jetIter];
      float trackMax=fForestJets.trackMax[jetIter];
      float trackSum=fForestJets.trackSum[jetIter];
      int trackN=fForestJets.trackN[jetIter];
      float trackHardSum=fForestJets.trackHardSum[jetIter];
      int trackHardN=fForestJets.trackHardN[jetIter];
      float chargedMax=fForestJets.chargedMax[jetIter];
      float chargedSum=fForestJets.chargedSum[jetIter];
      int chargedN=fForestJets.chargedN[jetIter];
      float chargedHardSum=fForestJets.chargedHardSum[jetIter];
      float chargedHardN=fForestJets.chargedHardN[jetIter];
      float photonMax=fForestJets.photonMax[jetIter];
      float photonSum=fForestJets.photonSum[jetIter];
      float photonHardSum=fForestJets.photonHardSum[jetIter];
      int photonHardN=fForestJets.photonHardN[jetIter];
      float neutralMax=fForestJets.neutralMax[jetIter];
      float neutralSum=fForestJets.neutralSum[jetIter];
      int neutralN=fForestJets.neutralN[jetIter];
      float hcalSum=fForestJets.hcalSum[jetIter];
      float ecalSum=fForestJets.ecalSum[jetIter];
      float eMax=fForestJets.eMax[jetIter];
      float eSum=fForestJets.eSum[jetIter];
      int eN=fForestJets.eN[jetIter];
      float muMax=fForestJets.muMax[jetIter];
      float muSum=fForestJets.muSum[jetIter];
      int muN=fForestJets.muN[jetIter];
      int   PfCHF=fForestJets.jtPfCHF[jetIter];
      float PfNHF=fForestJets.jtPfNHF[jetIter];
      float PfCEF=fForestJets.jtPfCEF[jetIter];
      float PfNEF=fForestJets.jtPfNEF[jetIter];
      int PfMUF=fForestJets.jtPfMUF[jetIter];
      int PfCHM=fForestJets.jtPfCHM[jetIter];
      int PfNHM=fForestJets.jtPfNHM[jetIter];
      int  PfCEM=fForestJets.jtPfCEM[jetIter];
      int PfNEM=fForestJets.jtPfNEM[jetIter];


      for(size_t ij=0; ij<tkJetsP4.size(); ij++) {
        if(jp4.DeltaR( tkJetsP4[ij] ) >0.4) continue;
        matchedJetsIdx.push_back(std::make_tuple(
						 ij,nsvtxTk,msvtx,csvVal,
						 phi, pu, mass,
						 trackMax,   trackSum,   trackN, trackHardSum, trackHardN,
						 chargedMax, chargedSum, chargedN, chargedHardSum, chargedHardN,
						 photonMax, photonSum, photonHardSum, photonHardN,
						 neutralMax,neutralSum,neutralN,
						 hcalSum,   ecalSum,
						 eMax,      eSum, eN, muMax,muSum, muN,
						 PfCHF, PfNHF, PfCEF, PfNEF, PfMUF, PfCHM, PfNHM, PfCEM, PfNEM
						 ));
	break;
      }


      if(jp4.Pt()<30.) continue;
      if(fabs(jp4.Eta())>2.4) continue;
      if(jp4.DeltaR(selLeptons[0])<0.4 || jp4.DeltaR(selLeptons[1])<0.4) continue;            
      bool isBTagged(csvVal>csvWP);
      
      pfJetsIdx.push_back(std::make_tuple(
					  jetIter,nsvtxTk,msvtx,csvVal,
					  phi, pu, mass,
					  trackMax,   trackSum,   trackN, trackHardSum, trackHardN,
					  chargedMax, chargedSum, chargedN, chargedHardSum, chargedHardN,
					  photonMax, photonSum, photonHardSum, photonHardN,
					  neutralMax,neutralSum,neutralN,
					  hcalSum,   ecalSum,
					  eMax,      eSum, eN, muMax,muSum, muN,
					  PfCHF, PfNHF, PfCEF, PfNEF, PfMUF, PfCHM, PfNHM, PfCEM, PfNEM
					  ));
      pfJetsP4.push_back(jp4);
      npfjets++;
      npfbjets += isBTagged;
      if(isBTagged && fabs(jp4.Eta())>1.2) allPFBInEB=false;

      if(npfjets==1){
        float dphi2ll(ll.DeltaPhi(jp4));
        if(fabs(dphi2ll)>2*TMath::Pi()/3.) hasAwayPFJet=true;
      }
    }
    std::sort(pfJetsIdx.begin(),      pfJetsIdx.end(),      orderByBtagInfo);


    //finalize analysing track jets
    std::sort(matchedJetsIdx.begin(), matchedJetsIdx.end(), orderByBtagInfo);
    bool allTkBInEB(true),hasAwayTkJet(false);
    float ntkjets(0),nbtkjets(0);
    for(size_t ij=0; ij<min(matchedJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<integral(JetInfo::index)>(matchedJetsIdx[ij]));
      float csv(std::get<integral(JetInfo::csvV2)>(matchedJetsIdx[ij]));
      TLorentzVector p4=tkJetsP4[idx];
      if(p4.Pt()<15) continue;
      bool isBTagged(csv>csvWP);
      ntkjets ++;
      nbtkjets += isBTagged;
      if(isBTagged && fabs(p4.Eta())>1.2) allTkBInEB=false;
      
      if(ntkjets==1){
        float dphi2ll(ll.DeltaPhi(p4));
        if(fabs(dphi2ll)>2*TMath::Pi()/3.) hasAwayTkJet=true;
      }
    }

    
    //fill control histograms
    std::vector<TString> categs;
    categs.push_back(dilCat);

    //monitor where the electrons are reconstructed
    if(dilCode==11*11 || dilCode==11*13){

      bool l1EE(fabs(selLeptons[0].Eta())>barrelEndcapEta[1]);
      bool l2EE(fabs(selLeptons[1].Eta())>barrelEndcapEta[1]);
      TString etaCateg(l2EE ? "E" : "B");
      if(dilCode==11*11) etaCateg += l1EE ? "E" : "B";

      categs.push_back(dilCat+etaCateg);
      if(isZ) {
        categs.push_back(dilCat+etaCateg+"Z");
        if(ntkjets==1 && hasAwayTkJet) {
          categs.push_back(dilCat+"Zawaytkj");
          categs.push_back(dilCat+etaCateg+"Zawaytkj");
        }
        if(npfjets==1 && hasAwayPFJet) {
          categs.push_back(dilCat+"Zawaypfj");
          categs.push_back(dilCat+etaCateg+"Zawaypfj");
        }
      }
    }
    

    std::vector<TString> addCategs;

    //monitor also after run where EE scale shift changed
    if(!isPP){
      addCategs.clear();
      TString pf( fForestTree.run>=firstEEScaleShiftRun ? "after" : "before" );
      for(auto c : categs) {
        addCategs.push_back(c); addCategs.push_back(c+pf); 
      }
      categs=addCategs;
    }

    //monitor according to the b-tagging category
    addCategs.clear();
    TString pfbcat(Form("%dpfb",min(npfbjets,2)));
    TString tkbcat(Form("%dpfb",min(npfbjets,2))); //fixme this should be based on the number of b-tagged track jets
    for(auto c : categs) { 
      addCategs.push_back(c); 
      addCategs.push_back(c+tkbcat); 
      addCategs.push_back(c+pfbcat); 
    }
    categs=addCategs;

    //monitor according to the centrality of the jets and electrons
    if(allEleInEB) {
      addCategs.clear();
      for(auto c: categs) {
        addCategs.push_back(c);
        if(allTkBInEB) addCategs.push_back(c+"alltkeb");
        if(allPFBInEB) addCategs.push_back(c+"allpfeb");
      }
    }


    float plotWgt(isMC ? fForestTree.weight : 1.0);
    for(int i=0; i<2; i++) {
      TString pf(Form("l%d",i+1));
      float pt(selLeptons[i].Pt());
      ht.fill(pf+"pt",        pt,                         plotWgt, categs);
      ht.fill(pf+"eta",       fabs(selLeptons[i].Eta()),  plotWgt, categs);
      float chiso(std::get<0>(liso[i]));
      float phoiso(std::get<1>(liso[i]));
      float neuiso(std::get<2>(liso[i]));
      ht.fill(pf+"chreliso",  chiso/pt,  plotWgt, categs);
      ht.fill(pf+"phoreliso", phoiso/pt,  plotWgt, categs);
      ht.fill(pf+"neureliso", neuiso/pt,  plotWgt, categs);
      ht.fill2D(pf+"chrelisovscen",  chiso/pt,   cenBin, plotWgt, categs);
      ht.fill2D(pf+"phorelisovscen", phoiso/pt,  cenBin, plotWgt, categs);
      ht.fill2D(pf+"neurelisovscen", neuiso/pt,  cenBin, plotWgt, categs);
    }

    ht.fill( "dphill",    fabs(selLeptons[0].DeltaPhi(selLeptons[1])), plotWgt, categs);
    ht.fill( "detall",    fabs(selLeptons[0].Eta()-selLeptons[1].Eta()), plotWgt, categs);
    ht.fill( "mll",       ll.M(),                                      plotWgt, categs);
    ht.fill( "ptll",      ll.Pt(),                                     plotWgt, categs);
    ht.fill( "npfjets",   npfjets,                                     plotWgt, categs);
    ht.fill( "npfbjets",  npfbjets,                                    plotWgt, categs);
    ht.fill( "ntkjets",   ntkjets,                                     plotWgt, categs);
    ht.fill( "ntkbjets",  nbtkjets,                                    plotWgt, categs);

    for(size_t ij=0; ij<min(matchedJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<integral(JetInfo::index)>(matchedJetsIdx[ij]));
      int ntks(std::get<integral(JetInfo::ntks)>(matchedJetsIdx[ij]));
      float svm(std::get<integral(JetInfo::msvtx)>(matchedJetsIdx[ij]));
      float csv(std::get<integral(JetInfo::csvV2)>(matchedJetsIdx[ij]));
      float phi(std::get<integral(JetInfo::phi)>(matchedJetsIdx[ij]));
      float pu(std::get<integral(JetInfo::pu)>(matchedJetsIdx[ij]));
      float mass(std::get<integral(JetInfo::mass)>(matchedJetsIdx[ij]));
      float trackMax(std::get<integral(JetInfo::trackMax)>(matchedJetsIdx[ij]));
      float trackSum(std::get<integral(JetInfo::trackSum)>(matchedJetsIdx[ij]));
      int trackN(std::get<integral(JetInfo::trackN)>(matchedJetsIdx[ij]));
      float trackHardSum(std::get<integral(JetInfo::trackHardSum)>(matchedJetsIdx[ij]));
      int trackHardN(std::get<integral(JetInfo::trackHardN)>(matchedJetsIdx[ij]));
      float chargedMax(std::get<integral(JetInfo::chargedMax)>(matchedJetsIdx[ij]));
      float chargedSum(std::get<integral(JetInfo::chargedSum)>(matchedJetsIdx[ij]));
      int chargedN(std::get<integral(JetInfo::chargedN)>(matchedJetsIdx[ij]));
      float chargedHardSum(std::get<integral(JetInfo::chargedHardSum)>(matchedJetsIdx[ij]));
      float chargedHardN(std::get<integral(JetInfo::chargedHardN)>(matchedJetsIdx[ij]));
      float photonMax(std::get<integral(JetInfo::photonMax)>(matchedJetsIdx[ij]));
      float photonSum(std::get<integral(JetInfo::photonSum)>(matchedJetsIdx[ij]));
      float photonHardSum(std::get<integral(JetInfo::photonHardSum)>(matchedJetsIdx[ij]));
      int photonHardN(std::get<integral(JetInfo::photonHardN)>(matchedJetsIdx[ij]));
      float neutralMax(std::get<integral(JetInfo::neutralMax)>(matchedJetsIdx[ij]));
      float neutralSum(std::get<integral(JetInfo::neutralSum)>(matchedJetsIdx[ij]));
      int neutralN(std::get<integral(JetInfo::neutralN)>(matchedJetsIdx[ij]));
      float hcalSum(std::get<integral(JetInfo::hcalSum)>(matchedJetsIdx[ij]));
      float ecalSum(std::get<integral(JetInfo::ecalSum)>(matchedJetsIdx[ij]));
      float eMax(std::get<integral(JetInfo::eMax)>(matchedJetsIdx[ij]));
      float eSum(std::get<integral(JetInfo::eSum)>(matchedJetsIdx[ij]));
      int eN(std::get<integral(JetInfo::eN)>(matchedJetsIdx[ij]));
      float muMax(std::get<integral(JetInfo::muMax)>(matchedJetsIdx[ij]));
      float muSum(std::get<integral(JetInfo::muSum)>(matchedJetsIdx[ij]));
      int muN(std::get<integral(JetInfo::muN)>(matchedJetsIdx[ij]));
      int   PfCHF(std::get<integral(JetInfo::PfCHF)>(matchedJetsIdx[ij]));
      float PfNHF(std::get<integral(JetInfo::PfNHF)>(matchedJetsIdx[ij]));
      float PfCEF(std::get<integral(JetInfo::PfCEF)>(matchedJetsIdx[ij]));
      float PfNEF(std::get<integral(JetInfo::PfNEF)>(matchedJetsIdx[ij]));
      int PfMUF(std::get<integral(JetInfo::PfMUF)>(matchedJetsIdx[ij]));
      int PfCHM(std::get<integral(JetInfo::PfCHM)>(matchedJetsIdx[ij]));
      int PfNHM(std::get<integral(JetInfo::PfNHM)>(matchedJetsIdx[ij]));
      int PfCEM(std::get<integral(JetInfo::PfCEM)>(matchedJetsIdx[ij]));
      int PfNEM(std::get<integral(JetInfo::PfNEM)>(matchedJetsIdx[ij]));

      TLorentzVector p4=tkJetsP4[idx];
      TString ppf(ij==1 ? "1" : "2");
      ht.fill( "tk"+ppf+"jbalance",  p4.Pt()/ll.Pt(),  plotWgt, categs);
      ht.fill( "tk"+ppf+"jpt",      p4.Pt(),          plotWgt, categs);
      ht.fill( "tk"+ppf+"jeta",     fabs(p4.Eta()),   plotWgt, categs);
      ht.fill( "tk"+ppf+"jsvtxm",   ntks,             plotWgt, categs);
      ht.fill( "tk"+ppf+"jsvtxntk", svm,              plotWgt, categs);
      ht.fill( "tk"+ppf+"jcsv",     csv,              plotWgt, categs);
      ht.fill( "tk"+ppf+"jphi",     phi,              plotWgt, categs);
      ht.fill( "tk"+ppf+"pu",       pu,               plotWgt, categs);
      ht.fill( "tk"+ppf+"mass",     mass,             plotWgt, categs);
      ht.fill( "tk"+ppf+"trackMax", trackMax,         plotWgt, categs);
      ht.fill( "tk"+ppf+"trackSum", trackSum,         plotWgt, categs);
      ht.fill( "tk"+ppf+"trackN",   trackN,           plotWgt, categs);
      ht.fill( "tk"+ppf+"trackHardSum",  trackHardSum,plotWgt, categs);
      ht.fill( "tk"+ppf+"trackHardN",    trackHardN, plotWgt, categs);
      ht.fill( "tk"+ppf+"chargedMax",    chargedMax, plotWgt, categs);
      ht.fill( "tk"+ppf+"chargedSum", chargedSum,     plotWgt, categs);
      ht.fill( "tk"+ppf+"chargedN", chargedN,         plotWgt, categs);
      ht.fill( "tk"+ppf+"chargedHardSum",  chargedHardSum,  plotWgt, categs);
      ht.fill( "tk"+ppf+"chargedHardN",   chargedHardN,     plotWgt, categs);
      ht.fill( "tk"+ppf+"photonMax",      photonMax,        plotWgt, categs);
      ht.fill( "tk"+ppf+"photonSum",     photonSum,         plotWgt, categs);
      ht.fill( "tk"+ppf+"photonHardSum", photonHardSum,     plotWgt, categs);
      ht.fill( "tk"+ppf+"photonHardN", photonHardN,         plotWgt, categs);
      ht.fill( "tk"+ppf+"neutralMax",   neutralMax,         plotWgt, categs);
      ht.fill( "tk"+ppf+"neutralSum",  neutralSum,          plotWgt, categs);
      ht.fill( "tk"+ppf+"neutralN",  neutralN,              plotWgt, categs);
      ht.fill( "tk"+ppf+"hcalSum",   hcalSum,               plotWgt, categs);
      ht.fill( "tk"+ppf+"ecalSum", ecalSum,                 plotWgt, categs);
      ht.fill( "tk"+ppf+"eMax", eMax,                       plotWgt, categs);
      ht.fill( "tk"+ppf+"eSum",  eSum,                      plotWgt, categs);
      ht.fill( "tk"+ppf+"eN", eN,                           plotWgt, categs);
      ht.fill( "tk"+ppf+"muMax", muMax,                     plotWgt, categs);
      ht.fill( "tk"+ppf+"muSum", muSum,                     plotWgt, categs);
      ht.fill( "tk"+ppf+"muN",   muN,                       plotWgt, categs);
      ht.fill( "tk"+ppf+"PfCHF",  PfCHF,                    plotWgt, categs);
      ht.fill( "tk"+ppf+"PfNHF",  PfNHF,                    plotWgt, categs);
      ht.fill( "tk"+ppf+"PfCEF",  PfCEF,                    plotWgt, categs);
      ht.fill( "tk"+ppf+"PfNEF", PfNEF,                     plotWgt, categs);
      ht.fill( "tk"+ppf+"PfMUF", PfMUF,                     plotWgt, categs);
      ht.fill( "tk"+ppf+"PfCHM",  PfCHM,                    plotWgt, categs);
      ht.fill( "tk"+ppf+"PfNHM", PfNHM,                     plotWgt, categs);
      ht.fill( "tk"+ppf+"PfCEM", PfCEM,                     plotWgt, categs);
      ht.fill( "tk"+ppf+"PfNEM", PfNEM,                     plotWgt, categs);

    }
    ht.fill( "tkrho", tkrho,            plotWgt, categs);
    
    for(size_t ij=0; ij<min(pfJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<integral(JetInfo::index)>(pfJetsIdx[ij]));
      int ntks(std::get<integral(JetInfo::ntks)>(pfJetsIdx[ij]));
      float svm(std::get<integral(JetInfo::msvtx)>(pfJetsIdx[ij]));
      float csv(std::get<integral(JetInfo::csvV2)>(pfJetsIdx[ij]));
      float phi(std::get<integral(JetInfo::phi)>(pfJetsIdx[ij]));
      float pu(std::get<integral(JetInfo::pu)>(pfJetsIdx[ij]));
      float mass(std::get<integral(JetInfo::mass)>(pfJetsIdx[ij]));
      float trackMax(std::get<integral(JetInfo::trackMax)>(pfJetsIdx[ij]));
      float trackSum(std::get<integral(JetInfo::trackSum)>(pfJetsIdx[ij]));
      int trackN(std::get<integral(JetInfo::trackN)>(pfJetsIdx[ij]));
      float trackHardSum(std::get<integral(JetInfo::trackHardSum)>(pfJetsIdx[ij]));
      int trackHardN(std::get<integral(JetInfo::trackHardN)>(pfJetsIdx[ij]));
      float chargedMax(std::get<integral(JetInfo::chargedMax)>(pfJetsIdx[ij]));
      float chargedSum(std::get<integral(JetInfo::chargedSum)>(pfJetsIdx[ij]));
      int chargedN(std::get<integral(JetInfo::chargedN)>(pfJetsIdx[ij]));
      float chargedHardSum(std::get<integral(JetInfo::chargedHardSum)>(pfJetsIdx[ij]));
      float chargedHardN(std::get<integral(JetInfo::chargedHardN)>(pfJetsIdx[ij]));
      float photonMax(std::get<integral(JetInfo::photonMax)>(pfJetsIdx[ij]));
      float photonSum(std::get<integral(JetInfo::photonSum)>(pfJetsIdx[ij]));
      float photonHardSum(std::get<integral(JetInfo::photonHardSum)>(pfJetsIdx[ij]));
      int photonHardN(std::get<integral(JetInfo::photonHardN)>(pfJetsIdx[ij]));
      float neutralMax(std::get<integral(JetInfo::neutralMax)>(pfJetsIdx[ij]));
      float neutralSum(std::get<integral(JetInfo::neutralSum)>(pfJetsIdx[ij]));
      int neutralN(std::get<integral(JetInfo::neutralN)>(pfJetsIdx[ij]));
      float hcalSum(std::get<integral(JetInfo::hcalSum)>(pfJetsIdx[ij]));
      float ecalSum(std::get<integral(JetInfo::ecalSum)>(pfJetsIdx[ij]));
      float eMax(std::get<integral(JetInfo::eMax)>(pfJetsIdx[ij]));
      float eSum(std::get<integral(JetInfo::eSum)>(pfJetsIdx[ij]));
      int eN(std::get<integral(JetInfo::eN)>(pfJetsIdx[ij]));
      float muMax(std::get<integral(JetInfo::muMax)>(pfJetsIdx[ij]));
      float muSum(std::get<integral(JetInfo::muSum)>(pfJetsIdx[ij]));
      int muN(std::get<integral(JetInfo::muN)>(pfJetsIdx[ij]));
      int   PfCHF(std::get<integral(JetInfo::PfCHF)>(pfJetsIdx[ij]));
      float PfNHF(std::get<integral(JetInfo::PfNHF)>(pfJetsIdx[ij]));
      float PfCEF(std::get<integral(JetInfo::PfCEF)>(pfJetsIdx[ij]));
      float PfNEF(std::get<integral(JetInfo::PfNEF)>(pfJetsIdx[ij]));
      int PfMUF(std::get<integral(JetInfo::PfMUF)>(pfJetsIdx[ij]));
      int PfCHM(std::get<integral(JetInfo::PfCHM)>(pfJetsIdx[ij]));
      int PfNHM(std::get<integral(JetInfo::PfNHM)>(pfJetsIdx[ij]));
      int PfCEM(std::get<integral(JetInfo::PfCEM)>(pfJetsIdx[ij]));
      int PfNEM(std::get<integral(JetInfo::PfNEM)>(pfJetsIdx[ij]));

      TLorentzVector p4=pfJetsP4[idx];
      TString ppf(ij==1 ? "1" : "2");
      ht.fill( "pf"+ppf+"jbalance",  p4.Pt()/ll.Pt(), plotWgt, categs);
      ht.fill( "pf"+ppf+"jpt",      p4.Pt(),         plotWgt, categs);
      ht.fill( "pf"+ppf+"jeta",     fabs(p4.Eta()),  plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxm",   ntks,            plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxntk", svm,             plotWgt, categs);
      ht.fill( "pf"+ppf+"jcsv",     csv,             plotWgt, categs);
      ht.fill( "pf"+ppf+"jphi",     phi,              plotWgt, categs);
      ht.fill( "pf"+ppf+"pu",       pu,               plotWgt, categs);
      ht.fill( "pf"+ppf+"mass",     mass,             plotWgt, categs);
      ht.fill( "pf"+ppf+"trackMax", trackMax,         plotWgt, categs);
      ht.fill( "pf"+ppf+"trackSum", trackSum,         plotWgt, categs);
      ht.fill( "pf"+ppf+"trackN",   trackN,           plotWgt, categs);
      ht.fill( "pf"+ppf+"trackHardSum",  trackHardSum,plotWgt, categs);
      ht.fill( "pf"+ppf+"trackHardN",    trackHardN, plotWgt, categs);
      ht.fill( "pf"+ppf+"chargedMax",    chargedMax, plotWgt, categs);
      ht.fill( "pf"+ppf+"chargedSum", chargedSum,     plotWgt, categs);
      ht.fill( "pf"+ppf+"chargedN", chargedN,         plotWgt, categs);
      ht.fill( "pf"+ppf+"chargedHardSum",  chargedHardSum,  plotWgt, categs);
      ht.fill( "pf"+ppf+"chargedHardN",   chargedHardN,     plotWgt, categs);
      ht.fill( "pf"+ppf+"photonMax",      photonMax,        plotWgt, categs);
      ht.fill( "pf"+ppf+"photonSum",     photonSum,         plotWgt, categs);
      ht.fill( "pf"+ppf+"photonHardSum", photonHardSum,     plotWgt, categs);
      ht.fill( "pf"+ppf+"photonHardN", photonHardN,         plotWgt, categs);
      ht.fill( "pf"+ppf+"neutralMax",   neutralMax,         plotWgt, categs);
      ht.fill( "pf"+ppf+"neutralSum",  neutralSum,          plotWgt, categs);
      ht.fill( "pf"+ppf+"neutralN",  neutralN,              plotWgt, categs);
      ht.fill( "pf"+ppf+"hcalSum",   hcalSum,               plotWgt, categs);
      ht.fill( "pf"+ppf+"ecalSum", ecalSum,                 plotWgt, categs);
      ht.fill( "pf"+ppf+"eMax", eMax,                       plotWgt, categs);
      ht.fill( "pf"+ppf+"eSum",  eSum,                      plotWgt, categs);
      ht.fill( "pf"+ppf+"eN", eN,                           plotWgt, categs);
      ht.fill( "pf"+ppf+"muMax", muMax,                     plotWgt, categs);
      ht.fill( "pf"+ppf+"muSum", muSum,                     plotWgt, categs);
      ht.fill( "pf"+ppf+"muN",   muN,                       plotWgt, categs);
      ht.fill( "pf"+ppf+"PfCHF",  PfCHF,                    plotWgt, categs);
      ht.fill( "pf"+ppf+"PfNHF",  PfNHF,                    plotWgt, categs);
      ht.fill( "pf"+ppf+"PfCEF",  PfCEF,                    plotWgt, categs);
      ht.fill( "pf"+ppf+"PfNEF", PfNEF,                     plotWgt, categs);
      ht.fill( "pf"+ppf+"PfMUF", PfMUF,                     plotWgt, categs);
      ht.fill( "pf"+ppf+"PfCHM",  PfCHM,                    plotWgt, categs);
      ht.fill( "pf"+ppf+"PfNHM", PfNHM,                     plotWgt, categs);
      ht.fill( "pf"+ppf+"PfCEM", PfCEM,                     plotWgt, categs);
      ht.fill( "pf"+ppf+"PfNEM", PfNEM,                     plotWgt, categs);

    }

  }

  //save histos to file  
  if(outURL!=""){
    TFile *fOut=TFile::Open(outURL,"RECREATE");
    fOut->cd();

    //store the weight sum for posterior normalization
    if(isMC) {
      TH1D *wgtH=new TH1D("wgtsum","wgtsum",1,0,1);
      wgtH->SetBinContent(1,wgtSum);
      wgtH->SetDirectory(fOut);
      wgtH->Write();
    }
    for (auto& it : ht.getPlots())  { 
      if(it.second->GetEntries()==0) continue;
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    for (auto& it : ht.get2dPlots())  { 
      if(it.second->GetEntries()==0) continue;
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    fOut->Close();
  }

  return 0;
}
