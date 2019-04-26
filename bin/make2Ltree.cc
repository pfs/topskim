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
#include "HeavyIonsAnalysis/topskim/include/ForestLeptons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestSkim.h"
#include "HeavyIonsAnalysis/topskim/include/ForestPFCands.h"
#include "HeavyIonsAnalysis/topskim/include/ForestJets.h"
#include "HeavyIonsAnalysis/topskim/include/LumiRun.h"
#include "HeavyIonsAnalysis/topskim/include/HistTool.h"
#include "HeavyIonsAnalysis/topskim/include/PFAnalysis.h"
#include "HeavyIonsAnalysis/topskim/include/LeptonSummary.h"
#include "HeavyIonsAnalysis/topskim/include/ForestGen.h"
#include "HeavyIonsAnalysis/topskim/include/ElectronId.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "../scripts/functions.cc"


const bool isDebug = true;
const float lepPtCut  = 20.;
const float lepEtaCut = 2.4;
//see https://indico.cern.ch/event/803679/contributions/3342407/attachments/1808912/2953435/egm-minipog-190308.pdf
const float eeScaleShift = 6.8182E-2/5.9097E-2;
const int firstEEScaleShiftRun = 327402; 
const float barrelEndcapEta[2]={1.4442,1.5660};
const float csvWP = 0.8838;

using namespace std;
using namespace fastjet;
using namespace TMVA;

//dispersion of rapidities of the final state
std::vector<float> getRapidityMoments(std::vector<TLorentzVector> & coll){
  std::vector<float> mom(3,0.);
  for(size_t i=0; i<coll.size(); i++) {
    TLorentzVector pi(coll[i]);
    mom[0]+=pi.Rapidity();
    mom[1]+=pow(pi.Rapidity(),2);
    for(size_t j=i; j<coll.size(); j++) {
      TLorentzVector pj(coll[j]);
      float dy=fabs(pj.Rapidity()-pi.Rapidity());
      mom[2]=max(dy,mom[2]);
    }
  }
  
  mom[0]=mom[0]/float(mom.size());
  mom[1]=sqrt(mom[1]/float(mom.size())-pow(mom[0],2));
  return mom;
}

int getRhoIndex(float eta){
        if      (eta < -3.0) return 1;
        else if (eta < -2.1) return 2;
        else if (eta < -1.3) return 2;
        else if (eta <  1.3) return 3;
        else if (eta <  2.1) return 4;
        else if (eta <  3.0) return 5;
        else return 6;
}


// index, ntks in svtx, m svtx, csv
typedef std::tuple<int,int,float,float,TLorentzVector,int,int> BtagInfo_t;
static bool orderByBtagInfo(const BtagInfo_t &a, const BtagInfo_t &b)
{
  //int ntks_a(std::get<1>(a)), ntks_b(std::get<1>(b));
  //if(ntks_a>ntks_b) return true;

  float csv_a(std::get<3>(a)), csv_b(std::get<3>(b));
  if(csv_a>csv_b) return true;
  return false;
}


static bool orderByPt(const LeptonSummary &a, const LeptonSummary &b)
{
  float pt_a(a.p4.Pt()), pt_b(b.p4.Pt());
  if(pt_a>pt_b) return true;
  return false;
}



//
int main(int argc, char* argv[])
{
  bool blind(true);
  TString inURL,outURL;
  bool isMC(false),isPP(false);
  int maxEvents(-1);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--in")!=string::npos && i+1<argc)       { inURL=TString(argv[i+1]); i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc) { outURL=TString(argv[i+1]); i++;}
    else if(arg.find("--max")!=string::npos && i+1<argc) { sscanf(argv[i+1],"%d",&maxEvents); }
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
  ht.addHist("fidcounter",  new TH2F("fidcounter", ";Fiducial counter;Events",5,0,5,1080,0,1080));
  ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(1,"all");
  ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(2,"=2l");
  ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(3,"=2l fid");
  ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(4,"=2l,#geq1b fid");
  ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(5,"=2l,#geq2b fid");


  if(!isMC) ht.addHist("ratevsrun",lumiTool.getLumiMonitor());

  //generic histograms
  for(int i=0; i<2; i++) {
    TString pf(Form("l%d",i+1));
    ht.addHist(pf+"pt",             new TH1F(pf+"pt",            ";Lepton transverse momentum [GeV];Events",20,20,200));
    ht.addHist(pf+"eta",            new TH1F(pf+"eta",           ";Lepton pseudo-rapidity;Events",20,0,2.5));

    for(int j=0; j<3; j++) {
      TString comp("ch");
      if(j==1) comp="pho";
      if(j==2) comp="nh";
      ht.addHist(pf+comp+"iso",      new TH1F(pf+comp+"iso",      ";PF "+comp+" isolation;Leptons",50,0,250));
      ht.addHist(pf+comp+"reliso",   new TH1F(pf+comp+"reliso",   ";Relative PF "+comp+" isolation;Leptons",50,0,2.0));
      ht.addHist(pf+comp+"isovscen", new TH2F(pf+comp+"isovscen", ";Centrality bin;PF "+comp+" isolation [GeV];Leptons",10,0,100,50,0,100));
      ht.addHist(pf+comp+"isovsrho", new TH2F(pf+comp+"isovsrho", ";#rho_{"+comp+"};PF "+comp+" isolation [GeV];Leptons",10,0,100,20,0,100));
    }    
  }

  //electron specific
  ht.addHist("esihih",      new TH1F("esihih",      ";#sigma(i#etai#eta);Electrons",       50,0,0.06));
  ht.addHist("edetavtx",    new TH1F("edetavtx",    ";#Delta#eta(vtx);Electrons",          50,0,0.015));
  ht.addHist("edphivtx",    new TH1F("edphivtx",    ";#Delta#phi(vtx) [rad];Electrons",    50,0,0.015));
  ht.addHist("ehoe",        new TH1F("ehoe"    ,    ";h/e;Electrons",                      50,0,0.25));
  ht.addHist("eempinv",     new TH1F("eempinv",     ";|1/E-1/p| [1/GeV];Electrons",        50,0,0.05));
  ht.addHist("ed0",         new TH1F("ed0",         ";d_{0} [cm];Electrons",               50,0,0.05));
  ht.addHist("edz",         new TH1F("edz",         ";d_{z} [cm];Electrons",               50,0,0.05));

  //muon specific
  ht.addHist("mmusta",     new TH1F("mmusta",      ";Muon stations;Muons",            15,0,15));   
  ht.addHist("mtrklay",    new TH1F("mtrklay",     ";Tracker layers;Muons",           25,0,25));
  ht.addHist("mchi2ndf",   new TH1F("mchi2ndf",    ";#chi^2/ndf;Muons",               50,0,15));
  ht.addHist("mmuhits",    new TH1F("mmuhits",     ";Muon hits;Muons",                25,0,25));
  ht.addHist("mpxhits",    new TH1F("mpxhits",     ";Pixel hits;Muons",               15,0,15));
  ht.addHist("md0",        new TH1F("md0",         ";d_{0} [cm];Muons",               50,0,0.5));
  ht.addHist("mdz",        new TH1F("mdz",         ";d_{z} [cm];Muons",               50,0,1.0));
 
  ht.addHist("mll",      new TH1F("mll",      ";Dilepton invariant mass [GeV];Events",40,0,200));
  ht.addHist("ptll",     new TH1F("ptll",     ";Dilepton transverse momentum [GeV];Events",25,0,200));
  ht.addHist("ptsum",    new TH1F("ptsum",    ";p_{T}(l)+p_{T}(l') [GeV];Events",25,0,200));
  ht.addHist("acopl",    new TH1F("acopl" ,   ";1-#Delta#phi(l,l')/#pi;Events",20,0,1.0));
  ht.addHist("detall",   new TH1F("detall",   ";#Delta#eta(l,l');Events",20,0,4));
  ht.addHist("drll",     new TH1F("drll"  ,   ";#DeltaR(l,l');Events",20,0,2*TMath::Pi()));

  ht.addHist("pfrapavg",      new TH1F("pfrapavg",     ";Average rapidity;Events",25,0,2.5));
  ht.addHist("pfraprms",      new TH1F("pfraprms",     ";#sigma(rapidity);Events",50,0,2.5));
  ht.addHist("pfrapmaxspan",  new TH1F("pfrapmaxspan", ";Maximum rapidity span;Events",25,0,5));
  ht.addHist("pfht",          new TH1F("pfht",         ";H_{T} [GeV];Events",25,0,500));
  ht.addHist("pfmht",         new TH1F("pfmht",        ";Missing H_{T} [GeV];Events",25,0,200));
  ht.addHist("npfjets",    new TH1F("npfjets",    ";Jet multiplicity;Events",8,0,8));
  ht.addHist("npfbjets",   new TH1F("npfbjets",   ";b-jet multiplicity;Events",5,0,5));    
  ht.addHist("npfsvtx",    new TH1F("npfsvtx",    ";Secondary vertex multiplicity;Events",5,0,5));
  for(size_t j=1; j<=2; j++){
    TString ppf(j==1 ? "1" : "2");
    ht.addHist("pf"+ppf+"jbalance",    new TH1F("pf"+ppf+"jbalance", ";R = p_{T}(j)/p_{T}(ll);Events",50,0,3));
    ht.addHist("pf"+ppf+"jpt",         new TH1F("pf"+ppf+"jpt",      ";Jet transverse momentum [GeV];Events",30,00,300));
    ht.addHist("pf"+ppf+"jeta",        new TH1F("pf"+ppf+"jeta",     ";Jet pseudo-rapidity;Events",20,0,2.5));
    ht.addHist("pf"+ppf+"jetavsphi",   new TH2F("pf"+ppf+"jetavsphi", ";Jet pseudo-rapidity;Jet azimuthal angle [rad];Events",100,-2.5,2.5,100,-TMath::Pi(),TMath::Pi()));
    ht.addHist("pf"+ppf+"jsvtxm",      new TH1F("pf"+ppf+"jsvtxm",   ";Secondary vertex mass;Events",25,0,6));
    ht.addHist("pf"+ppf+"jsvtxntk",    new TH1F("pf"+ppf+"jsvtxntk", ";Secondary vertex track multiplicity;Events",5,0,5));
    ht.addHist("pf"+ppf+"jcsv",        new TH1F("pf"+ppf+"jcsv",     ";CSVv2;Events",25,0,1));
  }

  //Get Tree info
  char* read = new char[100];
  TChain *hiInfoTree_p  = new TChain("HiForest/HiForestInfo");
  hiInfoTree_p->Add(inURL);
  TBranch *branch = hiInfoTree_p->GetBranch("GlobalTag");
  branch->SetAddress((void*)read);
  hiInfoTree_p->GetEntry(0);
  string GT(read, 0, 100);

  //Get global event filters 
  TChain *globalTree_p     = new TChain("skimanalysis/HltTree");
  globalTree_p->Add(inURL);
  ForestSkim fForestSkim(globalTree_p);

  //configure leptons
  TString lepTreeName("ggHiNtuplizerGED/EventTree");
  if(isPP) lepTreeName="ggHiNtuplizer/EventTree";
  if(GT.find("75X_mcRun2")!=string::npos) lepTreeName="ggHiNtuplizer/EventTree";
  TChain *lepTree_p     = new TChain(lepTreeName);
  lepTree_p->Add(inURL);
  ForestLeptons fForestLep(lepTree_p);
  ForestGen fForestGen(lepTree_p);

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
    TString muTrig("HLT_HIL3Mu20_v1");
    if( !hltTree_p->FindBranch(muTrig) ) muTrig="HLT_HIL3Mu15ForPPRef_v1";
    if( GT.find("75X_mcRun2")!=string::npos ) muTrig="HLT_HIL2Mu15ForPPRef_v1";
    hltTree_p->SetBranchStatus(muTrig,1);
    hltTree_p->SetBranchAddress(muTrig,&mtrig);
    TString eTrig("HLT_HIEle20_WPLoose_Gsf_v1");
    if( !hltTree_p->FindBranch(eTrig) ) eTrig="HLT_HISinglePhoton20_Eta3p1ForPPRef_v1";
    if( GT.find("75X_mcRun2")!=string::npos ) eTrig="HLT_HISinglePhoton40_Eta3p1ForPPRef_v1";
    hltTree_p->SetBranchStatus(eTrig,1);
    hltTree_p->SetBranchAddress(eTrig,&etrig);
    cout << "Using " << muTrig << " " << eTrig << " as MC triggers" << endl;
  }else{
    //hltTree_p->SetBranchStatus("HLT_HIL3Mu15_v1",1);
    //hltTree_p->SetBranchAddress("HLT_HIL3Mu15_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIL3Mu12_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu12_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20Gsf_v1",&etrig);    
  }

  TChain *rhoTree_p = new TChain("hiFJRhoAnalyzer/t");
  rhoTree_p->Add(inURL);
  std::vector<Double_t> *t_rho=0,*t_rhom=0;
  if(rhoTree_p){
    rhoTree_p->SetBranchAddress("rho", &t_rho);
    rhoTree_p->SetBranchAddress("rhom", &t_rhom);
  }else{
    std::cout << "[WARN] Can't find rho tree hiFJRhoAnalyzer/t" << std::endl;
  }

  
  // =============================================================
  // marc here make the output tree

  TTree * outTree = new TTree("tree", "tree with 2lepton selection and combined collections");

  // event and trigger variables
  Int_t  t_run, t_lumi, t_etrig, t_mtrig, t_isData;
  Long_t t_event;
  Float_t t_weight, t_cenbin;
  outTree->Branch("run"   , &t_run  , "run/I");
  outTree->Branch("lumi"  , &t_lumi , "lumi/I");
  outTree->Branch("event" , &t_event, "event/L");
  outTree->Branch("isData", &t_isData, "isData/I");

  outTree->Branch("weight", &t_weight, "weight/F");

  // centrality and different flavors of rho
  outTree->Branch("cenbin", &t_cenbin, "cenbin/F");
  
  outTree->Branch("rho",    &t_rho);
  outTree->Branch("rhom",   &t_rhom);

  Float_t t_globalrho;
  outTree->Branch("globalrho", &t_globalrho, "globalrho/F");

  outTree->Branch("etrig" , &t_etrig , "etrig/I");
  outTree->Branch("mtrig" , &t_mtrig , "mtrig/I");

  // variables per lepton, including iso
  Int_t t_nlep, t_lep_ind1, t_lep_ind2;
  std::vector<Float_t> t_lep_pt, t_lep_eta, t_lep_phi, t_lep_d0, t_lep_dz, t_lep_d0err, t_lep_phiso, t_lep_chiso, t_lep_nhiso, t_lep_rho, t_lep_isofull, t_lep_miniiso,t_lep_isofull20,t_lep_isofull25,t_lep_isofull30;
  std::vector<Bool_t> t_lep_matched;
  std::vector<Int_t  > t_lep_pdgId, t_lep_charge,t_lep_idflags;
  outTree->Branch("nlep"       , &t_nlep      , "nlep/I"            );
  outTree->Branch("lep_ind1"   , &t_lep_ind1  , "lep_ind1/I");
  outTree->Branch("lep_ind2"   , &t_lep_ind2  , "lep_ind2/I");
  outTree->Branch("lep_pt"     , &t_lep_pt);
  outTree->Branch("lep_eta"    , &t_lep_eta);
  outTree->Branch("lep_phi"    , &t_lep_phi);
  outTree->Branch("lep_idflags", &t_lep_idflags);
  outTree->Branch("lep_d0"    , &t_lep_d0);
  outTree->Branch("lep_d0err"  , &t_lep_d0err);
  outTree->Branch("lep_dz"     , &t_lep_dz);
  outTree->Branch("lep_phiso"  , &t_lep_phiso);
  outTree->Branch("lep_chiso"  , &t_lep_chiso);
  outTree->Branch("lep_nhiso"  , &t_lep_nhiso);
  outTree->Branch("lep_rho"    , &t_lep_rho);
  outTree->Branch("lep_pdgId"  , &t_lep_pdgId);
  outTree->Branch("lep_charge" , &t_lep_charge);
  outTree->Branch("lep_isofull", &t_lep_isofull);
  outTree->Branch("lep_isofull20", &t_lep_isofull20);
  outTree->Branch("lep_isofull25", &t_lep_isofull25);
  outTree->Branch("lep_isofull30", &t_lep_isofull30);
  outTree->Branch("lep_miniiso", &t_lep_miniiso);
  outTree->Branch("lep_matched", &t_lep_matched);

  // variables from dilepton system
  Float_t t_llpt, t_lleta, t_llphi, t_llm, t_dphi, t_deta, t_sumeta;
  outTree->Branch("llpt"   , &t_llpt   , "llpt/F");
  outTree->Branch("lleta"  , &t_lleta  , "lleta/F");
  outTree->Branch("llphi"  , &t_llphi  , "llphi/F");
  outTree->Branch("llm"    , &t_llm    , "llm/F");
  outTree->Branch("dphi"   , &t_dphi   , "dphi/F");
  outTree->Branch("deta"   , &t_deta   , "deta/F");
  outTree->Branch("sumeta" , &t_sumeta , "sumeta/F");

  // variables per jet (jets ordered by pt)
  Int_t t_njet;
  std::vector<Float_t> t_jet_pt, t_jet_eta, t_jet_phi, t_jet_mass, t_jet_csvv2;
  outTree->Branch("njet"      , &t_njet      , "njet/I"            );
  outTree->Branch("jet_pt"    , &t_jet_pt    , "jet_pt[njet]/F"    );
  outTree->Branch("jet_eta"   , &t_jet_eta   , "jet_eta[njet]/F"   );
  outTree->Branch("jet_phi"   , &t_jet_phi   , "jet_phi[njet]/F"   );
  outTree->Branch("jet_mass"  , &t_jet_mass  , "jet_mass[njet]/F" );
  outTree->Branch("jet_csvv2" , &t_jet_csvv2 , "jet_csvv2[njet]/F");

  // variables per bjet (jets ordered by csvv2)
  Int_t t_nbjet;
  std::vector<Float_t> t_bjet_pt, t_bjet_eta, t_bjet_phi, t_bjet_mass, t_bjet_csvv2;
  std::vector<Bool_t> t_bjet_drSafe;
  outTree->Branch("nbjet"      , &t_nbjet      , "nbjet/I"            );
  outTree->Branch("bjet_pt"    , &t_bjet_pt    );
  outTree->Branch("bjet_eta"   , &t_bjet_eta   );
  outTree->Branch("bjet_phi"   , &t_bjet_phi   );
  outTree->Branch("bjet_mass"  , &t_bjet_mass  );
  outTree->Branch("bjet_csvv2" , &t_bjet_csvv2 );
  outTree->Branch("bjet_drSafe" , &t_bjet_drSafe );

  std::vector<Float_t> t_bjet_matchpt, t_bjet_matcheta, t_bjet_matchphi, t_bjet_matchmass;
  outTree->Branch("bjet_genpt"    , &t_bjet_matchpt    );
  outTree->Branch("bjet_geneta"   , &t_bjet_matcheta   );
  outTree->Branch("bjet_genphi"   , &t_bjet_matchphi   );
  outTree->Branch("bjet_genmass"  , &t_bjet_matchmass  );

  std::vector<Int_t> t_bjet_flavor, t_bjet_flavorForB;
  outTree->Branch("bjet_flavor"  , &t_bjet_flavor  );
  outTree->Branch("bjet_flavorB"  , &t_bjet_flavorForB  );

  // constructed variables like ht and stuff
  Float_t t_ht, t_mht, t_apt, t_dphilll2;
  outTree->Branch("ht"     , &t_ht     , "ht/F");
  outTree->Branch("mht"    , &t_mht    , "mht/F");
  outTree->Branch("apt"    , &t_apt    , "apt/F");
  outTree->Branch("dphilll2"    , &t_dphilll2    , "dphilll2/F");

  Float_t t_bdt, t_bdt_rarity, t_fisher2;
  outTree->Branch("bdt", &t_bdt, "bdt/F");
  outTree->Branch("bdtrarity", &t_bdt_rarity, "bdtrarity/F");
  outTree->Branch("fisher2", &t_fisher2, "fisher2/F");
  // =============================================================
  //
  TMVA::Tools::Instance();
  TMVA::Reader *reader        = new TMVA::Reader( "!Color:!Silent" );
  TMVA::Reader *readerFisher2 = new TMVA::Reader( "!Color:!Silent" );

  // make new variables because i'm too lazy to think right now
  Float_t bdt_l1pt, bdt_apt, bdt_abslleta, bdt_dphilll2, bdt_sumabseta;//, bdt_flavor;

  // these must have the same name as in the training. and the same ordeeeeeeer
  // copy directly from the script that runs the training:
  //dataloader.AddVariable('lep_pt[0]'  , 'p_{T}^{lep1}'     , 'GeV' , 'F')
  //dataloader.AddVariable('apt'        , 'A_{pt}'           , ''    , 'F')
  //dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'       , 'GeV' , 'F')
  //dataloader.AddVariable('abs(lleta)' , '|#eta^{ll}|'      , ''    , 'F')
  //dataloader.AddVariable('dphi'       , '|#Delta #phi|'    , 'rad' , 'F')
  //dataloader.AddVariable('abs(lep_eta[0])+abs(lep_eta[1])' , '#sum |#eta_{i}|', ''    , 'F')
  //dataloader.AddVariable('abs(lep_pdgId[0]*lep_pdgId[1])'  , 'flavor', ''    , 'F')
  //
  reader->AddVariable("lep_pt[0]"  , &bdt_l1pt    );
  reader->AddVariable("apt"        , &bdt_apt     );
  reader->AddVariable("llpt"       , &t_llpt      );
  reader->AddVariable("abs(lleta)" , &bdt_abslleta);
  reader->AddVariable("abs(dphi)"  , &t_dphi      );
  reader->AddVariable("abs(lep_eta[0])+abs(lep_eta[1])", &bdt_sumabseta);

  // for the fisher just take these two
  //dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'       , 'GeV' , 'F')
  //dataloader.AddVariable('dphi'       , '|#Delta #phi|'    , 'rad' , 'F')
  readerFisher2->AddVariable("llpt", &t_llpt);
  readerFisher2->AddVariable("abs(dphi)", &t_dphi);

  TString methodName       ("BDTG method");
  TString methodNameFisher2("Fisher method");
  // hard coded path for now ...
  TString weightFile("/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/HeavyIonsAnalysis/topskim/scripts/training_dy/weights/TMVAClassification_BDTG.weights.xml");
  reader->BookMVA( methodName, weightFile);

  TString weightFileFisher2("/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/HeavyIonsAnalysis/topskim/scripts/training_dy_fisher2/weights/TMVAClassification_Fisher.weights.xml");
  readerFisher2->BookMVA( methodNameFisher2, weightFileFisher2);

    
  Double_t wgtSum(0);
  std::vector<Double_t> allWgtSum;
  int nEntries = (int)lepTree_p->GetEntries();  
  int entryDiv = ((int)(nEntries/20));    
  cout << inURL << " has " << nEntries << " events to process" << endl;
  if(maxEvents>0) { 
    nEntries=TMath::Min(nEntries,maxEvents); 
    cout << "Number of events to process limited to " << nEntries << endl;
  }
  for(int entry = 0; entry < nEntries; entry++){
    
    if(entry%entryDiv == 0) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;

    globalTree_p->GetEntry(entry);
    lepTree_p->GetEntry(entry);
    pfCandTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);    
    hltTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    if(rhoTree_p) rhoTree_p->GetEntry(entry);

    //gen level analysis
    float evWgt(1.0);
    int genDileptonCat(1.);
    std::vector<TLorentzVector> genLeptons, genBjets;
    bool isGenDilepton(false),isLeptonFiducial(false),is1bFiducial(false),is2bFiducial(false);    
    if(isMC) {
     
      //gen level selection      
      for(size_t i=0; i<fForestGen.mcPID->size(); i++) {
        int pid=fForestGen.mcPID->at(i);
        //int sta=fForestGen.mcStatus->at(i);
        int mom_pid=fForestGen.mcMomPID->at(i);
        int gmom_pid=fForestGen.mcGMomPID->at(i);
        
        if( abs(pid)<6  && abs(mom_pid)==6 ) {
          TLorentzVector p4(0,0,0,0);
          p4.SetPtEtaPhiM( fForestGen.mcPt->at(i), fForestGen.mcEta->at(i), fForestGen.mcPhi->at(i), fForestGen.mcMass->at(i) );
          if(p4.Pt()>30 && fabs(p4.Eta())<2.5) genBjets.push_back(p4);          
        }

        //leptons from t->W->l or W->tau->l
        if( abs(pid)==11 || abs(pid)==13 ) {
          
          bool isFromW( abs(mom_pid)==24 && abs(gmom_pid)==6 );
          bool isTauFeedDown( abs(mom_pid)==15 && abs(gmom_pid)==24 );
          if(isFromW || isTauFeedDown) {
            TLorentzVector p4(0,0,0,0);
            p4.SetPtEtaPhiM( fForestGen.mcPt->at(i), fForestGen.mcEta->at(i), fForestGen.mcPhi->at(i), fForestGen.mcMass->at(i) );
            //if(p4.Pt()>20 && fabs(p4.Eta())<2.5) 
            genLeptons.push_back(p4);
            genDileptonCat *= abs(pid);
          }
        }
      }
      
      isGenDilepton=(genLeptons.size()==2);      
      isLeptonFiducial=(isGenDilepton && 
                        genLeptons[0].Pt()>20 && fabs(genLeptons[0].Eta())<2.5 && 
                        genLeptons[1].Pt()>20 && fabs(genLeptons[1].Eta())<2.5);
      is1bFiducial=(isLeptonFiducial && genBjets.size()>0);
      is2bFiducial=(isLeptonFiducial && genBjets.size()>1);
      
      //event weights and fiducial counters   
      if(isMC) {
        evWgt=fForestTree.ttbar_w->size()==0 ? 1. : fForestTree.ttbar_w->at(0);
        size_t nWgts(fForestTree.ttbar_w->size());
        if(nWgts==0) nWgts=1;
        if(allWgtSum.size()==0) allWgtSum.resize(nWgts,0.);
        for(size_t i=0; i<nWgts; i++) {
          Double_t iwgt(fForestTree.ttbar_w->size()==0 ? 1. : fForestTree.ttbar_w->at(i));
          allWgtSum[i]+=iwgt;
          ht.fill2D("fidcounter",0,i,iwgt,"gen");
          if(isGenDilepton)    ht.fill2D("fidcounter",1,i,iwgt,"gen");
          if(isLeptonFiducial) ht.fill2D("fidcounter",2,i,iwgt,"gen");
          if(is1bFiducial)     ht.fill2D("fidcounter",3,i,iwgt,"gen");
          if(is2bFiducial)     ht.fill2D("fidcounter",4,i,iwgt,"gen");
        }
      }
    }
    
    wgtSum += evWgt;    
    float plotWgt(evWgt);
    
    //first of all require a trigger
    int trig=etrig+mtrig;
    if(trig==0) continue;
    if(isSingleMuPD) {
      if(mtrig==0) continue;
      if(etrig!=0) continue;
    }
    if(isSingleElePD) {
      if(etrig==0) continue;
    }

    //apply global filters
    if(!isMC && GT.find("103X")!=string::npos){
      if(TMath::Abs(fForestTree.vz) > 20) continue;
      if(!fForestSkim.phfCoincFilter2Th4) continue;
      if(!fForestSkim.pclusterCompatibilityFilter) continue;
      if(!fForestSkim.pprimaryVertexFilter) continue;
    }
    else if(!isMC && GT.find("75X")!=string::npos){
      if(TMath::Abs(fForestTree.vz) > 15) continue;
      if(!fForestSkim.phfCoincFilter) continue;
      if(!fForestSkim.HBHENoiseFilterResult) continue;
      if(!fForestSkim.pcollisionEventSelection) continue;
      if(!fForestSkim.pprimaryVertexFilter) continue;
    }

    //build jets from different PF candidate collections   
    SlimmedPFCollection_t pfColl;
    for(size_t ipf=0; ipf<fForestPF.pfId->size(); ipf++) {
      int id(abs(fForestPF.pfId->at(ipf)));
      float mass(0.13957);  //pions
      if(id==4) mass=0.;    //photons
      if(id>=5) mass=0.497; //K0L
      pfColl.push_back( getSlimmedPF( id, fForestPF.pfPt->at(ipf),fForestPF.pfEta->at(ipf),fForestPF.pfPhi->at(ipf),mass) );
    }

    Float_t globalrho = getRho(pfColl,{1,2,3,4,5,6},-1.,5.);

    //monitor trigger and centrality
    float cenBin=0;
    bool isCentralEvent(false);
    if(!isMC){
      cenBin=0.5*fForestTree.hiBin;
      isCentralEvent=(!isPP && cenBin<30);
      Int_t runBin=lumiTool.getRunBin(fForestTree.run);
      Float_t lumi=lumiTool.getLumi(fForestTree.run);
      if(lumi>0.){
        if(etrig>0) ht.fill("ratevsrun",runBin,1./lumi,"e");
        if(mtrig>0) ht.fill("ratevsrun",runBin,1./lumi,"m");
      }
    }

    //the selected leptons
    std::vector<LeptonSummary> selLeptons;

    //select muons
    std::vector<LeptonSummary> noIdMu;
    for(unsigned int muIter = 0; muIter < fForestLep.muPt->size(); ++muIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestLep.muPt->at(muIter),fForestLep.muEta->at(muIter),fForestLep.muPhi->at(muIter),0.1057);
      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(p4.Pt() < lepPtCut) continue;

      LeptonSummary l(13,p4);
      l.charge  = fForestLep.muCharge->at(muIter);
      l.chiso   = fForestLep.muPFChIso->at(muIter);
      l.nhiso   = fForestLep.muPFNeuIso->at(muIter);
      l.phoiso  = fForestLep.muPFPhoIso->at(muIter);
      int   tmp_rhoind  = getRhoIndex(p4.Eta());
      if (!isMC && GT.find("75X")==string::npos){
        float tmp_rho_par = 0.0013 * TMath::Power(t_rho->at(tmp_rhoind)+15.83,2) + 0.29 * (t_rho->at(tmp_rhoind)+15.83); 
        l.isofull = (l.chiso+l.nhiso+l.phoiso - tmp_rho_par)/p4.Pt();
      }
      else {
        l.isofull = -1.;
	  }
      l.rho = isPP ? globalrho : t_rho->at(tmp_rhoind);

      l.isofullR=getIsolationFull( pfColl, l.p4);
      l.miniiso = getMiniIsolation( pfColl ,l.p4, l.id);
      l.d0      = fForestLep.muD0   ->at(muIter);
      l.d0err   = 0.; //fForestLep.muD0Err->at(muIter); // no d0err for muons!!!
      l.dz      = fForestLep.muDz   ->at(muIter);
      l.origIdx = muIter;
      l.isMatched=false;
      for(size_t ig=0;ig<genLeptons.size(); ig++) {
        if(genLeptons[ig].DeltaR(l.p4)<0.1) continue;
        l.isMatched=true;
      }

      noIdMu.push_back(l);

      //id (Tight muon requirements)
      int type=fForestLep.muType->at(muIter);
      bool isGlobal( ((type>>1)&0x1) );
      if(!isGlobal) continue;
      bool isPF( ((type>>5)&0x1) );
      if(!isPF) continue;
      bool isGlobalMuonPromptTight(fForestLep.muChi2NDF->at(muIter)<10. && fForestLep.muMuonHits->at(muIter)>0);
      if(!isGlobalMuonPromptTight) continue;
      if(fForestLep.muStations->at(muIter)<=1) continue;
      if(fForestLep.muTrkLayers->at(muIter) <= 5) continue;
      if(fForestLep.muPixelHits->at(muIter) == 0) continue;
      if(TMath::Abs(fForestLep.muInnerD0->at(muIter)) >=0.2 ) continue;
      if(TMath::Abs(fForestLep.muInnerDz->at(muIter)) >=0.5) continue;

      l.idFlags=1;

      //selected a good muon
      selLeptons.push_back(l);
    }
    std::sort(noIdMu.begin(),noIdMu.end(),orderByPt);


    //monitor muon id variables
    if(noIdMu.size()>1) {
      TLorentzVector p4[2] = {noIdMu[0].p4,      noIdMu[1].p4};
      int midx[2]          = {noIdMu[0].origIdx, noIdMu[1].origIdx};
      int charge(noIdMu[0].charge*noIdMu[1].charge);
      if( fabs((p4[0]+p4[1]).M()-91)<15 ) {
        TString cat("zmmctrl");
        if(charge>0) cat="ss"+cat;
        for(size_t i=0; i<2; i++) {
          ht.fill("mmusta",    fForestLep.muStations->at(midx[i]),            plotWgt,cat);
          ht.fill("mtrklay",   fForestLep.muTrkLayers->at(midx[i]),           plotWgt,cat);
          ht.fill("mchi2ndf",  fForestLep.muChi2NDF->at(midx[i]),             plotWgt,cat);
          ht.fill("mmuhits",   fForestLep.muMuonHits->at(midx[i]),            plotWgt,cat);
          ht.fill("mpxhits",   fForestLep.muPixelHits->at(midx[i]),           plotWgt,cat);
          ht.fill("md0",       TMath::Abs(fForestLep.muInnerD0->at(midx[i])), plotWgt,cat);
          ht.fill("mdz",       TMath::Abs(fForestLep.muInnerDz->at(midx[i])), plotWgt,cat);          
        }
      }
    }

    //select electrons
    //cf. https://twiki.cern.ch/twiki/pub/CMS/HiHighPt2019/HIN_electrons2018_followUp.pdf
    std::vector<LeptonSummary> noIdEle;
    for(unsigned int eleIter = 0; eleIter < fForestLep.elePt->size(); ++eleIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestLep.elePt->at(eleIter),fForestLep.eleEta->at(eleIter),fForestLep.elePhi->at(eleIter),0.000511);

      //apply ad-hoc shift for endcap electrons if needed, i.e., PromptReco'18
      if(!isMC && fForestTree.run<=firstEEScaleShiftRun && TMath::Abs(p4.Eta())>=barrelEndcapEta[1] && GT.find("fixEcalADCToGeV")==string::npos && GT.find("75X")==string::npos)
        p4 *=eeScaleShift;         


      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(TMath::Abs(p4.Eta()) > barrelEndcapEta[0] && TMath::Abs(p4.Eta()) < barrelEndcapEta[1] ) continue;
      if(p4.Pt() < lepPtCut) continue;	      
      
      LeptonSummary l(11,p4);
      l.charge  = fForestLep.eleCharge->at(eleIter);
      if(GT.find("75X_mcRun2")==string::npos) {
	l.chiso   = fForestLep.elePFChIso03->at(eleIter);
	l.nhiso   = fForestLep.elePFNeuIso03->at(eleIter);
	l.phoiso  = fForestLep.elePFPhoIso03->at(eleIter);
      }else{
        l.chiso   = fForestLep.elePFChIso->at(eleIter);
        l.nhiso   = fForestLep.elePFNeuIso->at(eleIter);
        l.phoiso  = fForestLep.elePFPhoIso->at(eleIter);
      }
      int   tmp_rhoind  = getRhoIndex(p4.Eta());
      if (!isMC && GT.find("75X")==string::npos){
        float tmp_rho_par = 0.0011 * TMath::Power(t_rho->at(tmp_rhoind)+142.4,2) - 0.14 * (t_rho->at(tmp_rhoind)+142.4); 
        l.isofull = (l.chiso+l.nhiso+l.phoiso - tmp_rho_par)/p4.Pt();
      }
      else {
        l.isofull = -1.;
	  }

      l.rho = isPP ? globalrho : t_rho->at(tmp_rhoind);

      l.isofullR= getIsolationFull( pfColl, l.p4);
      l.miniiso = getMiniIsolation( pfColl ,l.p4, l.id);
      l.d0      = fForestLep.eleD0   ->at(eleIter);
      l.d0err   = fForestLep.eleD0Err->at(eleIter);
      l.dz      = fForestLep.eleDz   ->at(eleIter);
      l.origIdx=eleIter;
      l.isMatched=false;
      for(size_t ig=0;ig<genLeptons.size(); ig++) {
        if(genLeptons[ig].DeltaR(l.p4)<0.1) continue;
        l.isMatched=true;
      }

      noIdEle.push_back(l);
      
      l.idFlags=getElectronId(TMath::Abs(fForestLep.eleSCEta->at(eleIter))< barrelEndcapEta[0],
                                 fForestLep.eleSigmaIEtaIEta->at(eleIter),
                                 fForestLep.eledEtaAtVtx->at(eleIter),
                                 fForestLep.eledPhiAtVtx->at(eleIter),
                                 fForestLep.eleHoverE->at(eleIter),
                                 fForestLep.eleEoverPInv->at(eleIter),
                                 fForestLep.eleD0->at(eleIter),
                                 fForestLep.eleDz->at(eleIter),
                                 fForestLep.eleMissHits->at(eleIter),
                                 isCentralEvent);
      
      //id'ed electron
      if(!isLooseElectron(l.idFlags)) continue;
      selLeptons.push_back(l);
    }
    std::sort(noIdEle.begin(),noIdEle.end(),orderByPt);       

    //monitor electron id variables
    if(noIdEle.size()>1) {
      TLorentzVector p4[2] = {noIdEle[0].p4,noIdEle[1].p4};
      int eidx[2]          = {noIdEle[0].origIdx,noIdEle[1].origIdx};
      int charge(noIdEle[0].charge*noIdEle[1].charge);
      if( fabs((p4[0]+p4[1]).M()-91)<15 ) {
        TString basecat("zeectrl");
        if(charge>0) basecat="ss"+basecat;
        for(size_t i=0; i<2; i++) {
          TString cat(basecat);
          cat += (fabs(p4[i].Eta())>=barrelEndcapEta[1] ? "EE" : "EB");
          ht.fill("esihih",  fForestLep.eleSigmaIEtaIEta->at(eidx[i]),         plotWgt,cat);
          ht.fill("edetavtx", TMath::Abs(fForestLep.eledEtaAtVtx->at(eidx[i])), plotWgt,cat);
          ht.fill("edphivtx", TMath::Abs(fForestLep.eledPhiAtVtx->at(eidx[i])), plotWgt,cat);
          ht.fill("ehoe",     fForestLep.eleHoverE->at(eidx[i]),                plotWgt,cat);
          ht.fill("eempinv",  fForestLep.eleEoverPInv->at(eidx[i]),             plotWgt,cat);
          ht.fill("ed0",      TMath::Abs(fForestLep.eleD0->at(eidx[i])),        plotWgt,cat);
          ht.fill("edz",      TMath::Abs(fForestLep.eleDz->at(eidx[i])),        plotWgt,cat);          
        }
      }
    }

    //sort selected electrons by pt
    std::sort(selLeptons.begin(),selLeptons.end(),orderByPt);

    //apply basic preselection
    if(mtrig+etrig==0) continue;
    if(selLeptons.size()<2) continue;
    TLorentzVector ll(selLeptons[0].p4+selLeptons[1].p4);
    t_llpt=ll.Pt();
    t_lleta=ll.Eta();
    t_llphi=ll.Phi();
    t_llm=ll.M();
    t_dphi=TMath::Abs(selLeptons[0].p4.DeltaPhi(selLeptons[1].p4));
    t_deta=fabs(selLeptons[0].p4.Eta()-selLeptons[1].p4.Eta());
    t_sumeta=selLeptons[0].p4.Eta()+selLeptons[1].p4.Eta();
    int dilCode(selLeptons[0].id*selLeptons[1].id);
    TString dilCat("mm");
    if(dilCode==11*13) dilCat="em";
    if(dilCode==11*11) dilCat="ee";

    //ee and mm events should come from the appropriate primary dataset
    if(!isMC) {
      if(dilCode==11*11 && !isSingleElePD) continue;
      if(dilCode==13*13 && !isSingleMuPD) continue;
    }

    if(blind) {
      bool isZ( dilCode!=11*13 && fabs(t_llm-91)<15);
      int charge(selLeptons[0].charge*selLeptons[1].charge);
      if(!isMC && !isZ && charge<0 && fForestTree.run>=326887) continue;
    }      
          
    //analyze jets
    std::vector<bool> drSafe_pfJet;
    std::vector<BtagInfo_t> pfJetsIdx,nodr_pfJetsIdx;
    std::vector<TLorentzVector> pfJetsP4,nodr_pfJetsP4,nodr_pfJetsP4GenMatch;
    int npfjets(0),npfbjets(0); 
    for(int jetIter = 0; jetIter < fForestJets.nref; jetIter++){

      //at least two tracks
      if(fForestJets.trackN[jetIter]<2) continue;

      TLorentzVector jp4(0,0,0,0);
      jp4.SetPtEtaPhiM( fForestJets.jtpt[jetIter],fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);

      float csvVal=fForestJets.discr_csvV2[jetIter];
      int nsvtxTk=fForestJets.svtxntrk[jetIter];
      float msvtx=fForestJets.svtxm[jetIter];

      if(jp4.Pt()<30.) continue;
      if(fabs(jp4.Eta())>2.4) continue;
      bool isBTagged(csvVal>csvWP);      

      // simple matching to the closest jet in dR. require at least dR < 0.3
      TLorentzVector matchjp4(0,0,0,0);
      int refFlavor(0),refFlavorForB(0);
      if (isMC){
        //std::vector<TLorentzVector> matchedJets;
        for (int genjetIter = 0; genjetIter < fForestJets.ngen; genjetIter++){
          if (jetIter == fForestJets.genmatchindex[genjetIter]) {
            matchjp4.SetPtEtaPhiM( fForestJets.genpt[genjetIter],fForestJets.geneta[genjetIter],fForestJets.genphi[genjetIter],fForestJets.genm[genjetIter]);
          }
        }  
        refFlavor=fForestJets.refparton_flavor[jetIter];
        refFlavorForB=fForestJets.refparton_flavorForB[jetIter];
      }

      nodr_pfJetsIdx.push_back( std::make_tuple(nodr_pfJetsP4.size(),nsvtxTk,msvtx,csvVal,matchjp4,refFlavor,refFlavorForB) );
      nodr_pfJetsP4.push_back(jp4);
      bool isdrSafe(false);

      //cross clean wrt to leptons
      if(jp4.DeltaR(selLeptons[0].p4)<0.4 || jp4.DeltaR(selLeptons[1].p4)<0.4) {
        pfJetsIdx.push_back(std::make_tuple(pfJetsP4.size(),nsvtxTk,msvtx,csvVal,matchjp4,refFlavor,refFlavorForB));
        pfJetsP4.push_back(jp4);
        npfjets++;
        npfbjets += isBTagged;
        isdrSafe=true;
      }

      drSafe_pfJet.push_back(isdrSafe);      
    }
    std::sort(pfJetsIdx.begin(),       pfJetsIdx.end(),      orderByBtagInfo);
    std::sort(nodr_pfJetsIdx.begin(),  nodr_pfJetsIdx.end(), orderByBtagInfo);


    //for gen fill again fiducial counters
    if(isMC) {      
      
      bool isMatchedDilepton(abs(genDileptonCat)==abs(dilCode));
      if( (genDileptonCat==11*11 && etrig==0) || (genDileptonCat==13*13 && mtrig==0)) 
        isMatchedDilepton=false;

      bool isIsoDilepton(true);
      for(size_t i=0; i<2; i++) {
        if( (abs(selLeptons[i].id)==13 && selLeptons[i].isofull<0.26) ||
            (abs(selLeptons[i].id)==11 && selLeptons[i].isofull<0.16) ) continue;
        isIsoDilepton=false;
      }

      std::vector<TString> fidCats;
      fidCats.push_back( isMatchedDilepton   ? "lep"    : "fakelep" );
      if(isIsoDilepton) {
        fidCats.push_back( isMatchedDilepton ? "isolep" : "fakeisolep" );
        if(npfbjets>0) {
          fidCats.push_back( isMatchedDilepton && is1bFiducial ? "isolep1b" : "fakeisolep1b" );
          if(npfbjets>1) {
            fidCats.push_back( isMatchedDilepton && is2bFiducial ? "isolep2b" : "fakeisolep2b" );
          }
        }
      }
      
      size_t nWgts(fForestTree.ttbar_w->size());
      if(nWgts==0) nWgts=1;
      for(size_t i=0; i<nWgts; i++){
        Double_t iwgt(fForestTree.ttbar_w->size()==0 ? 1. : fForestTree.ttbar_w->at(i));
        ht.fill2D("fidcounter",0,i,iwgt,fidCats);
        if(isGenDilepton)    ht.fill2D("fidcounter",1,i,iwgt,fidCats);
        if(isLeptonFiducial) ht.fill2D("fidcounter",2,i,iwgt,fidCats);
        if(is1bFiducial)     ht.fill2D("fidcounter",3,i,iwgt,fidCats);
        if(is2bFiducial)     ht.fill2D("fidcounter",4,i,iwgt,fidCats);
      }
    }


    //define categories for pre-selection control histograms
    std::vector<TString> categs;
    categs.push_back(dilCat);
    
    std::vector<TString> addCategs;
    
    //monitor after run where EE scale shift changed
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
    for(auto c : categs) { 
      addCategs.push_back(c); 
      if(npfbjets==0) addCategs.push_back(c+"0pfb"); 
      if(npfbjets>0) addCategs.push_back(c+"geq1pfb"); 
    }
    categs=addCategs;

    
    //fill histogams
    for(int i=0; i<2; i++) {
      TString pf(Form("l%d",i+1));
      float pt(selLeptons[i].p4.Pt());
      ht.fill(pf+"pt",             pt,                            plotWgt, categs);
      ht.fill(pf+"eta",            fabs(selLeptons[i].p4.Eta()),  plotWgt, categs);
      ht.fill(pf+"chiso",          selLeptons[i].chiso,           plotWgt, categs);
      ht.fill(pf+"phoiso",         selLeptons[i].phoiso,          plotWgt, categs);
      ht.fill(pf+"nhiso",          selLeptons[i].nhiso,           plotWgt, categs);

    }

    ht.fill( "acopl",     1-fabs(t_dphi)/TMath::Pi(),                   plotWgt, categs);
    ht.fill( "detall",    t_deta,                                       plotWgt, categs);
    ht.fill( "drll",      selLeptons[0].p4.DeltaR(selLeptons[1].p4),    plotWgt, categs);
    ht.fill( "mll",       t_llm,                                        plotWgt, categs);
    ht.fill( "ptll",      t_llpt,                                       plotWgt, categs);
    ht.fill( "ptsum",     selLeptons[0].p4.Pt()+selLeptons[1].p4.Pt(),  plotWgt, categs);

    //PF jets
    ht.fill( "npfjets",   npfjets,   plotWgt, categs);
    ht.fill( "npfbjets",  npfbjets,  plotWgt, categs);
    std::vector<TLorentzVector> pfFinalState;
    pfFinalState.push_back(selLeptons[0].p4);
    pfFinalState.push_back(selLeptons[1].p4);
    for(size_t ij=0; ij<min(pfJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<0>(pfJetsIdx[ij]));
      int ntks(std::get<1>(pfJetsIdx[ij]));
      float svm(std::get<2>(pfJetsIdx[ij]));
      float csv(std::get<3>(pfJetsIdx[ij]));
      TLorentzVector p4=pfJetsP4[idx];
      if(csv>csvWP) pfFinalState.push_back(p4);
      TString ppf(ij==0 ? "1" : "2");
      ht.fill( "pf"+ppf+"jbalance", p4.Pt()/ll.Pt(), plotWgt, categs);
      ht.fill( "pf"+ppf+"jpt",      p4.Pt(),         plotWgt, categs);
      ht.fill( "pf"+ppf+"jeta",     fabs(p4.Eta()),  plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxm",   ntks,            plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxntk", svm,             plotWgt, categs);
      ht.fill( "pf"+ppf+"jcsv",     csv,             plotWgt, categs);
      ht.fill2D( "pf"+ppf+"jetavsphi",   p4.Eta(),p4.Phi(),   plotWgt, categs);
    }

    std::vector<float> rapMoments=getRapidityMoments(pfFinalState);
    ht.fill( "pfrapavg",     rapMoments[0], plotWgt, categs);
    ht.fill( "pfraprms",     rapMoments[1], plotWgt, categs);
    ht.fill( "pfrapmaxspan", rapMoments[2], plotWgt, categs);
    float pfht(0.);
    TLorentzVector vis(0,0,0,0);
    for(auto p : pfFinalState) { vis+=p; pfht+=p.Pt(); }
    ht.fill( "pfht",         pfht, plotWgt, categs);
    ht.fill( "pfmht",        vis.Pt(), plotWgt, categs);

    // for tree filling set all the proper variables
    t_run    = fForestTree.run;
    t_lumi   = fForestTree.lumi;
    t_event  = fForestTree.evt;
    t_weight = plotWgt;
    t_cenbin = cenBin;
    t_globalrho = globalrho;
    t_etrig  = etrig;
    t_mtrig  = mtrig;

    // fill the leptons ordered by pt
    t_lep_pt    .clear();
    t_lep_eta   .clear();
    t_lep_phi   .clear();
    t_lep_pdgId .clear();
    t_lep_idflags.clear();
    t_lep_d0 .clear();
    t_lep_d0err .clear();
    t_lep_dz  .clear();
    t_lep_charge.clear();
    t_lep_chiso.clear();    
    t_lep_phiso.clear();
    t_lep_nhiso.clear();
    t_lep_rho.clear();    
    t_lep_isofull.clear();
    t_lep_isofull20.clear();
    t_lep_isofull25.clear();
    t_lep_isofull30.clear();
    t_lep_miniiso.clear();
    t_lep_matched.clear();
    t_nlep = selLeptons.size();
    t_lep_ind1 = -1;
    t_lep_ind2 = -1;
    for (int ilep = 0; ilep < t_nlep; ++ilep){
      t_lep_pt    .push_back( selLeptons[ilep].p4.Pt()  );
      t_lep_eta   .push_back( selLeptons[ilep].p4.Eta() );
      t_lep_phi   .push_back( selLeptons[ilep].p4.Phi() );
      t_lep_idflags.push_back(selLeptons[ilep].idFlags);
      t_lep_d0    .push_back( selLeptons[ilep].d0 );
      t_lep_d0err .push_back( selLeptons[ilep].d0err);
      t_lep_dz    .push_back( selLeptons[ilep].dz  );
      t_lep_chiso .push_back( selLeptons[ilep].chiso );
      t_lep_phiso .push_back( selLeptons[ilep].phoiso );
      t_lep_nhiso .push_back( selLeptons[ilep].nhiso );
      t_lep_rho   .push_back( selLeptons[ilep].rho );
      t_lep_pdgId .push_back( selLeptons[ilep].id );
      t_lep_charge.push_back( selLeptons[ilep].charge );
      t_lep_isofull.push_back( selLeptons[ilep].isofull );
      t_lep_isofull20.push_back( selLeptons[ilep].isofullR[0] );
      t_lep_isofull25.push_back( selLeptons[ilep].isofullR[1] );
      t_lep_isofull30.push_back( selLeptons[ilep].isofullR[2] );
      t_lep_miniiso.push_back( selLeptons[ilep].miniiso );
      t_lep_matched.push_back( selLeptons[ilep].isMatched );
      if(selLeptons[ilep].isofull < 0.16 && t_lep_ind1 < 0) t_lep_ind1 = ilep;
      if(selLeptons[ilep].isofull < 0.16 && t_lep_ind1 > -0.5 && t_lep_ind2 < 0) t_lep_ind2 = ilep;
    }

    // fill the jets ordered by b-tag
    t_bjet_pt   .clear();
    t_bjet_eta  .clear();
    t_bjet_phi  .clear();
    t_bjet_mass .clear();
    t_bjet_csvv2.clear();
    t_bjet_drSafe.clear();
    t_bjet_matchpt  .clear();
    t_bjet_matcheta .clear();
    t_bjet_matchphi .clear();
    t_bjet_matchmass.clear();
    t_bjet_flavor.clear();
    t_bjet_flavorForB.clear();
    t_nbjet = nodr_pfJetsIdx.size();
    for (int ij = 0; ij < t_nbjet; ij++) {
      int idx = std::get<0>(nodr_pfJetsIdx[ij]);
      t_bjet_pt   .push_back( nodr_pfJetsP4[idx].Pt()  );
      t_bjet_eta  .push_back( nodr_pfJetsP4[idx].Eta() );
      t_bjet_phi  .push_back( nodr_pfJetsP4[idx].Phi() );
      t_bjet_mass .push_back( nodr_pfJetsP4[idx].M()   );
      t_bjet_csvv2.push_back( std::get<3>(nodr_pfJetsIdx[ij])   );
      t_bjet_drSafe.push_back( drSafe_pfJet[idx] );
      t_bjet_matchpt  .push_back( std::get<4>(nodr_pfJetsIdx[ij]).Pt());
      t_bjet_matcheta .push_back( std::get<4>(nodr_pfJetsIdx[ij]).Eta());
      t_bjet_matchphi .push_back( std::get<4>(nodr_pfJetsIdx[ij]).Phi());
      t_bjet_matchmass.push_back( std::get<4>(nodr_pfJetsIdx[ij]).M());
      t_bjet_flavor.push_back( std::get<5>(nodr_pfJetsIdx[ij]) );
      t_bjet_flavorForB.push_back( std::get<6>(nodr_pfJetsIdx[ij]) );
    }


    t_ht  = pfht;
    t_mht = vis.Pt();

    // now set the 4 variables that we added for the tmva reader for the bdt evaluation
    bdt_l1pt      = t_lep_pt[0];
    bdt_apt       = (t_lep_pt[0]-t_lep_pt[1])/(t_lep_pt[0]+t_lep_pt[1]);
    bdt_abslleta  = fabs(t_lleta);
    bdt_dphilll2  = fabs(dphi_2(t_lep_pt[0],t_lep_eta[0],t_lep_phi[0],t_lep_pt[1],t_lep_eta[1],t_lep_phi[1],2)); // this function is in functions.cc in scripts/
    bdt_sumabseta = fabs(t_lep_eta[0])+fabs(t_lep_eta[1]);
    //bdt_flavor    = abs(t_lep_pdgId[0]*t_lep_pdgId[1]); //abs should be fine here, it's an int
    t_apt         = bdt_apt;
    t_dphilll2    = bdt_dphilll2;
    t_bdt         = reader->EvaluateMVA( methodName );
    t_bdt_rarity  = reader->GetRarity  ( methodName );
    t_fisher2     = readerFisher2->EvaluateMVA( methodNameFisher2 );

    t_isData = !isMC;

    outTree->Fill();
  }

  //save histos to file  
  if(outURL!=""){
    TFile *fOut=TFile::Open(outURL,"RECREATE");
    fOut->cd();

    outTree->Write();

    //store the weight sum for posterior normalization
    TH1D *wgtH=new TH1D("wgtsum","wgtsum",1,0,1);
    wgtH->SetBinContent(1,wgtSum);
    wgtH->SetDirectory(fOut);
    wgtH->Write();

    TH1D *allwgtH=new TH1D("allwgtsum","allwgtsum",allWgtSum.size(),0,allWgtSum.size());
    for(size_t i=0; i<allWgtSum.size(); i++)
      allwgtH->SetBinContent(i+1,allWgtSum[i]);
    allwgtH->SetDirectory(fOut);
    allwgtH->Write();
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
