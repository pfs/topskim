#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TSystem.h"
#include "TF1.h"
#include "TRandom3.h"

#include <string>
#include <vector>

#include "HeavyIonsAnalysis/topskim/include/tnp_weight.h"
#include "HeavyIonsAnalysis/topskim/include/tnp_electrons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestHiTree.h"
#include "HeavyIonsAnalysis/topskim/include/ForestHLTObject.h"
#include "HeavyIonsAnalysis/topskim/include/ForestLeptons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestSkim.h"
#include "HeavyIonsAnalysis/topskim/include/ForestPFCands.h"
#include "HeavyIonsAnalysis/topskim/include/ForestJets.h"
#include "HeavyIonsAnalysis/topskim/include/JetCorrector.h"
#include "HeavyIonsAnalysis/topskim/include/JetUncertainty.h"
#include "HeavyIonsAnalysis/topskim/include/LumiRun.h"
#include "HeavyIonsAnalysis/topskim/include/HistTool.h"
#include "HeavyIonsAnalysis/topskim/include/PFAnalysis.h"
#include "HeavyIonsAnalysis/topskim/include/LeptonSummary.h"
#include "HeavyIonsAnalysis/topskim/include/ForestGen.h"
#include "HeavyIonsAnalysis/topskim/include/ElectronId.h"
#include "HeavyIonsAnalysis/topskim/include/BtagUncertaintyComputer.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "../scripts/functions.cc"


const bool isDebug = true;
const float lepPtCut  = 15.;
const float muEtaCut = 2.4;
const float eleEtaCut = 2.1;
//see https://indico.cern.ch/event/803679/contributions/3342407/attachments/1808912/2953435/egm-minipog-190308.pdf
const float eeScaleShift = 6.8182E-2/5.9097E-2;
const int firstEEScaleShiftRun = 327402; 
const float barrelEndcapEta[2]={1.4442,1.5660};
const float hem1516Eta[2]={-3.0,-1.9};
const float hem1516Phi[2]={-1.6,-0.9};
const float csvWPList[2]={0.81,0.91};
int csvWP = 0;
// runs with HLT issues: any path using tracking was disabled => includes L3 muon paths. 
// total lumi in these runs is 30.103/ub (total lumi available in golden json unblinde period is 425.349/ub)
std::vector<int> badMuonTriggerRuns={326482,326483,326500,326520,326527,326528,326530,326532,326533,326534,326535,326546,326548,326549,326550,326568,326569,326571};

using namespace std;
using namespace fastjet;

TRandom3 *smearRand = new TRandom3(42);

//Madgraph BR(W->lnu) correction
float getMadgraphBRWlCorrection(int nl) {
  if(nl==0) return 1.0236380625;
  if(nl==1) return 0.987973875;
  if(nl==2) return 0.95355225;
  return 1.0;
}

//eta-based indexing for rho
int getRhoIndex(float eta,std::vector<Double_t> *etaMin=NULL, std::vector<Double_t> *etaMax=NULL){
  if(etaMin==NULL || etaMax==NULL) return -1;
  for(size_t i=0; i<etaMin->size(); i++) {
    if(eta>=etaMin->at(i) && eta<etaMax->at(i)) return i;
  }
  return -1;
}

//defines a relativistic breit-wigner
TF1 *getRBW(float m,float g) {
  TF1 *bwigner=new TF1("bwigner",
                       "[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))",
                       max(float(0.),m-50*g),m+50*g);

  bwigner->SetParName(0,"N");
  bwigner->SetParameter(0,1.0);
  bwigner->SetParName(1,"m_{0}");
  bwigner->FixParameter(1,m);
  bwigner->SetParName(2,"#Gamma_{t}");
  bwigner->FixParameter(2,g);

  return bwigner;
}

//gets the weight for different (mass,width) hypotheses
float weightBW(TF1 *bwigner,std::vector<float> obsm,float g,float m,float gini,float mini) {

  //
  if(bwigner==NULL || obsm.size()!=2 || obsm[0]<150 || obsm[1]<150) return 1.;

  bwigner->FixParameter(1,mini);
  bwigner->FixParameter(2,gini);
  float nini=bwigner->Integral(max(m-50*g,float(0.)),m+50*g);
      
  bwigner->FixParameter(1,m);
  bwigner->FixParameter(2,g);
  float n=bwigner->Integral(max(m-50*g,float(0.)),m+50*g);

  float wgt(1.0);
  for(auto obsm_i : obsm){
    bwigner->FixParameter(1,mini);
    bwigner->FixParameter(2,gini);
    float vini=bwigner->Eval(obsm_i);

    bwigner->FixParameter(1,m);
    bwigner->FixParameter(2,g);
    float v=bwigner->Eval(obsm_i);

    wgt *= (v/n) / (vini/nini);
  }

  return wgt;
}


//electron calibration
float calibratedPt(float pt, float eta, float cen, bool isMC) {

  float newpt;
  float scale = 0.;

  if (isMC){
    if (fabs(eta) < 1.45) {
      if      (cen < 10.) scale = 0.974;
      else if (cen < 30.) scale = 0.992;
      else                scale = 1.005;
    }else {
      if      (cen < 10.) scale = 0.913;
      else if (cen < 30.) scale = 0.952;
      else                scale = 0.992;
    }
  } else {
    if (fabs(eta) < 1.45) {
      if      (cen < 10.) scale = 0.990;
      else if (cen < 30.) scale = 1.006;
      else                scale = 1.016;
    }else {
      if      (cen < 10.) scale = 0.976;
      else if (cen < 30.) scale = 1.015;
      else                scale = 1.052;
    }
  }

  newpt = pt * scale;

  float smear = 1.;

  if (isMC){
    if (fabs(eta) < 1.45) {
      if      (cen < 10.) smear = 0.904;
      else if (cen < 30.) smear = 1.379;
      else                smear = 1.786;
    } else {
      if      (cen < 10.) smear = 3.051;
      else if (cen < 30.) smear = 1.214;
      else                smear = 3.451;
    }

    newpt = newpt * smearRand->Gaus(1., smear / 91.1876);
  }

  return newpt;
}

// ordering based on b-tag info (index, ntks in svtx, m svtx, csv)
typedef std::tuple<int,int,float,float,TLorentzVector,int,int> BtagInfo_t;
static bool orderByBtagInfo(const BtagInfo_t &a, const BtagInfo_t &b)
{
  float csv_a(std::get<3>(a)), csv_b(std::get<3>(b));
  if(csv_a>csv_b) return true;
  return false;
}

//ordering based on pt
static bool orderByPt(const LeptonSummary &a, const LeptonSummary &b)
{
  float pt_a(a.p4.Pt()), pt_b(b.p4.Pt());
  if(pt_a>pt_b) return true;
  return false;
}

// btag efficiencies from the AN
float btagEfficiencies(int flavor, float cenbin){
  if(csvWP==0){
    if      (fabs(flavor) == 5) return (cenbin <= 30 ? 0.556 : 0.683); // bs
    else if (fabs(flavor) == 0) return (cenbin <= 30 ? 0.057 : 0.042); // unmatched
    else                        return (cenbin <= 30 ? 0.023 : 0.008); // udsg
  }
  else{
    if      (fabs(flavor) == 5) return (cenbin <= 30 ? 0.385 : 0.546); // bs
    else if (fabs(flavor) == 0) return (cenbin <= 30 ? 0.013 : 0.008); // unmatched
    else                        return (cenbin <= 30 ? 0.005 : 0.001); // udsg
  }
}

//parameterizations from Glauber MC
Double_t findNcoll(int hiBin) {
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
};


//
int main(int argc, char* argv[])
{
  bool blind(true);
  TString inURL,outURL;
  bool isMC(false),isPP(false),isAMCATNLO(false);
  int maxEvents(-1);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--in")!=string::npos && i+1<argc)         { inURL=TString(argv[i+1]); i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc)   { outURL=TString(argv[i+1]); i++;}
    else if(arg.find("--max")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&maxEvents); }
    else if(arg.find("--csvWP")!=string::npos && i+1<argc) { sscanf(argv[i+1],"%d",&csvWP); }
    else if(arg.find("--mc")!=string::npos)                { isMC=true;  }
    else if(arg.find("--pp")!=string::npos)                { isPP=true;  }
    else if(arg.find("--amcatnlo")!=string::npos)          { isAMCATNLO=true;  }
  }
  
  bool isSingleMuPD( !isMC && inURL.Contains("SkimMuons"));
  bool isSingleElePD( !isMC && inURL.Contains("SkimElectrons"));
  bool isMuSkimedMCPD( isMC && inURL.Contains("HINPbPbAutumn18DR_skims") && inURL.Contains("Muons"));
  bool isEleSkimedMCPD( isMC && inURL.Contains("HINPbPbAutumn18DR_skims") && inURL.Contains("Electrons"));
  LumiRun lumiTool;
  ElectronEfficiencyWrapper eleEff("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data", false);

  //read expected trigger efficiencies
  TString trigEffURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/trigeff_mc.root");
  gSystem->ExpandPathName(trigEffURL);
  TFile *fIn=TFile::Open(trigEffURL);
  TGraphAsymmErrors *e_mctrigeff=(TGraphAsymmErrors *)fIn->Get("e_pt_trigeff");
  fIn->Close();

  //read isolation efficiencies
  TString isosfURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/isolation_sf.root");
  gSystem->ExpandPathName(isosfURL);
  fIn=TFile::Open(isosfURL);
  std::map<TString, TH2 *> isoEffSFs;
  TString isomaps[]={"cen_169","periph_169","cen_121","periph_121"};
  for(size_t i=0; i<sizeof(isomaps)/sizeof(TString); i++){
    TString key(isomaps[i]);
    isoEffSFs[key]=(TH2*)fIn->Get("sfiso2eff_"+key);
    isoEffSFs[key]->SetDirectory(0);
  }
  fIn->Close();
  
  // initialize the JEC and associated unc files
  std::vector<std::string> FilesData;
  TString DATA_L2RelativeURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V4_DATA_L2Relative_AK4PF.txt");
  gSystem->ExpandPathName(DATA_L2RelativeURL);
  TString DATA_L2ResidualURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V4_DATA_L2Residual_AK4PF.txt");
  gSystem->ExpandPathName(DATA_L2ResidualURL);
  FilesData.push_back(DATA_L2RelativeURL.Data());
  FilesData.push_back(DATA_L2ResidualURL.Data());
  
  JetCorrector JECData(FilesData);
  TString JEUDataURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V4_DATA_Uncertainty_AK4PF.txt");
  gSystem->ExpandPathName(DATA_L2RelativeURL);
  JetUncertainty JEUData(DATA_L2RelativeURL.Data());
  
  std::vector<std::string> FilesMC;
  TString FilesMCURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V4_MC_L2Relative_AK4PF.txt");
  gSystem->ExpandPathName(FilesMCURL);
  
  FilesMC.push_back(FilesMCURL.Data());
  
  JetCorrector JECMC(FilesMC);
  TString JECMCURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V4_MC_Uncertainty_AK4PF.txt");
  gSystem->ExpandPathName(JECMCURL);
  JetUncertainty JEUMC(JECMCURL.Data());

  //define the quenching model parameterization
  TF1 * quenchingModel = new TF1("quenchingModel", "[0]/(TMath::Sqrt(2.*TMath::Pi())*0.73*x)*TMath::Exp(-1.*TMath::Power(TMath::Log(x/[0])+1.5,2)/2./0.73/0.73)", 0., 50.);
  TF1 * centralityModel = new TF1("centralityModel", "gaus");
  centralityModel->SetParameter(0,  1.090);
  centralityModel->SetParameter(1, -0.144);
  centralityModel->SetParameter(2,  0.442);

  if(isPP)
    cout << "Treating as a pp collision file" << endl;
  if(isMC)
    cout << "Treating as a MC file" << endl;
  if(isAMCATNLO)
    cout << "This is an amc@NLO file (for the ME weights)" << endl;
  else
    cout << "Will assume that this is a powheg file!" << endl;

  //weights of interest to store
  std::vector<size_t> meIdxList={0,1,2,3,4,6,8}; //qcd weights
  if(!isMC) 
    meIdxList.clear();
  else {
    if(isPP){
      for(size_t i=10; i<=111; i++) meIdxList.push_back(i); //hessian NNPDF3.1
      meIdxList.push_back(116); meIdxList.push_back(117);   //alphaS variation +/-0.001
    }else if(isAMCATNLO){    
      meIdxList[0]=1080;
      for(size_t i=1081; i<=1176; i++) meIdxList.push_back(i); //EPPS16nlo_CT14nlo_Pb208
      for(size_t i=1177; i<=1209; i++) meIdxList.push_back(i); //nCTEQ15FullNuc_208_82
    }else{    
      meIdxList[0]=864;
      for(size_t i=865; i<=960; i++) meIdxList.push_back(i); //EPPS16nlo_CT14nlo_Pb208
      for(size_t i=961; i<=993; i++) meIdxList.push_back(i); //nCTEQ15FullNuc_208_82
    }
    cout << "Will store " <<  meIdxList.size() << " ME weights" << endl;
  }

  //for breit-wigner reweithing
  TF1 *rbwigner( isMC ? getRBW(172.5,1.31) : NULL );

  // Initialize the btagging SF stuff
  BTagSFUtil * myBTagUtil = new BTagSFUtil(42);

  //attach t the relevant trees
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
  TChain *jetTree_p     = new TChain(isPP ? "ak4PFJetAnalyzer/t" : "akFlowPuCs4PFJetAnalyzer/t");
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
  TString muTrigName(""),eTrigName("");
  if(isPP){
    muTrigName="HLT_HIL3Mu20_v";
    if( !hltTree_p->FindBranch(muTrigName+"1") ) muTrigName="HLT_HIL3Mu15ForPPRef_v";
    if( GT.find("75X_mcRun2")!=string::npos ) muTrigName="HLT_HIL2Mu15ForPPRef_v";
    hltTree_p->SetBranchStatus(muTrigName+"1",1);
    hltTree_p->SetBranchAddress(muTrigName+"1",&mtrig);
    eTrigName="HLT_HIEle20_WPLoose_Gsf_v";
    if( !hltTree_p->FindBranch(eTrigName+"1") ) eTrigName="HLT_HISinglePhoton20_Eta3p1ForPPRef_v";
    if( GT.find("75X_mcRun2")!=string::npos ) eTrigName="HLT_HISinglePhoton40_Eta3p1ForPPRef_v";
    hltTree_p->SetBranchStatus(eTrigName+"1",1);
    hltTree_p->SetBranchAddress(eTrigName+"1",&etrig);
  }else{
    muTrigName="HLT_HIL3Mu12_v";
    eTrigName="HLT_HIEle20Gsf_v";
    if(isSingleMuPD || isMuSkimedMCPD) mtrig = 1;
    if(isSingleElePD || isEleSkimedMCPD) etrig = 1;
    if(isMC and !isMuSkimedMCPD and !isEleSkimedMCPD){
      hltTree_p->SetBranchStatus(muTrigName+"1",1);
      hltTree_p->SetBranchAddress(muTrigName+"1",&mtrig);
      hltTree_p->SetBranchStatus(eTrigName+"1",1);
      hltTree_p->SetBranchAddress(eTrigName+"1",&etrig);
    }
  }

  //trigger objects
  cout << "Using " << muTrigName << " " << eTrigName << " as triggers" << endl;
  TChain *muHLTObj_p =NULL, *eleHLTObj_p=NULL;
  ForestHLTObject *muHLTObjs=NULL, *eleHLTObjs=NULL;
  if(!isPP){
    muHLTObj_p=new TChain("hltobject/"+muTrigName);
    muHLTObj_p->Add(inURL);
    muHLTObjs=new ForestHLTObject(muHLTObj_p);
    eleHLTObj_p = new TChain("hltobject/"+eTrigName);
    eleHLTObj_p->Add(inURL);
    eleHLTObj_p->AddFriend(muHLTObj_p);
    eleHLTObjs=new ForestHLTObject(eleHLTObj_p);
  }

  TChain *rhoTree_p = new TChain("hiFJRhoAnalyzerFinerBins/t");
  rhoTree_p->Add(inURL);
  std::vector<Double_t> *t_rho=0,*t_rhom=0,*t_etaMin=0,*t_etaMax=0;
  if(rhoTree_p){
    rhoTree_p->SetBranchAddress("rho", &t_rho);
    rhoTree_p->SetBranchAddress("etaMin", &t_etaMin);
    rhoTree_p->SetBranchAddress("etaMax", &t_etaMax);
    rhoTree_p->SetBranchAddress("rhom", &t_rhom);
  }else{
    std::cout << "[WARN] Can't find rho tree hiFJRhoAnalyzerFinerBins/t" << std::endl;
  }
  
  
  //prepare output tree
  TTree * outTree = new TTree("tree", "tree with single lepton + jets selection");

  //global event variables
  Int_t  t_run, t_lumi, t_isData;
  Long_t t_event;
  Float_t t_vx, t_vy, t_vz, t_trigSF, t_trigSFUnc;
  std::vector<Float_t> t_meWeights;
  Float_t t_cenbin, t_ncollWgt;
  outTree->Branch("run"   , &t_run  ,  "run/I");
  outTree->Branch("lumi"  , &t_lumi ,  "lumi/I");
  outTree->Branch("event" , &t_event,  "event/L");
  outTree->Branch("isData", &t_isData, "isData/I");
  outTree->Branch("vx",     &t_vx,     "vx/F");
  outTree->Branch("vy",     &t_vy,     "vy/F");
  outTree->Branch("vz",     &t_vz,     "vz/F");
  outTree->Branch("rho",    &t_rho);
  outTree->Branch("rhom",   &t_rhom);
  outTree->Branch("cenbin",   &t_cenbin,  "cenbin/F");
  outTree->Branch("ncollWgt", &t_ncollWgt, "ncollWgt/F");

  //MC truth
  Bool_t t_genFiducial;
  Float_t t_weight, t_weight_BRW;
  Int_t t_genlep_id;
  Float_t t_genlep_pt,t_genlep_eta,t_genlep_phi;
  std::vector<Float_t> t_genjet_pt,t_genjet_eta,t_genjet_phi,t_genjet_m;
  outTree->Branch("weight",    &t_weight,     "weight/F");
  outTree->Branch("weightBRW", &t_weight_BRW, "weightBRW/F");
  outTree->Branch("meWeights", &t_meWeights);
  outTree->Branch("genFiducial" ,   &t_genFiducial ,"genFiducial/O");
  outTree->Branch("genlep_id" ,     &t_genlep_id ,  "genlep_id/I");
  outTree->Branch("genlep_pt" ,     &t_genlep_pt ,  "genlep_pt/F");
  outTree->Branch("genlep_eta" ,    &t_genlep_eta , "genlep_eta/F");
  outTree->Branch("genlep_phi" ,    &t_genlep_phi , "genlep_phi/F");
  outTree->Branch("genjet_pt" ,     &t_genjet_pt);
  outTree->Branch("genjet_eta" ,    &t_genjet_eta);
  outTree->Branch("genjet_phi" ,    &t_genjet_phi);
  outTree->Branch("genjet_m" ,      &t_genjet_m);

  //pf analysis
  Float_t t_globalrho;
  outTree->Branch("globalrho", &t_globalrho, "globalrho/F");

  //lepton variables
  Float_t t_lep_pt, t_lep_eta, t_lep_phi, t_lep_d0, t_lep_dz, t_lep_d0err, t_lep_phiso, t_lep_chiso, t_lep_nhiso, t_lep_rho, t_lep_isofull, t_lep_miniiso,t_lep_isofull20,t_lep_isofull25,t_lep_isofull30;
  Int_t   t_lep_pdgId, t_lep_charge,t_lep_idflags;
  Float_t t_lepSF,t_lepSFUnc,t_lepIsoSF,t_lepIsoSFUnc;  
  outTree->Branch("lep_pt"     ,   &t_lep_pt);
  outTree->Branch("lep_eta"    ,   &t_lep_eta);
  outTree->Branch("lep_phi"    ,   &t_lep_phi);
  outTree->Branch("lep_idflags",   &t_lep_idflags);
  outTree->Branch("lep_d0"     ,   &t_lep_d0);
  outTree->Branch("lep_d0err"  ,   &t_lep_d0err);
  outTree->Branch("lep_dz"     ,   &t_lep_dz);
  outTree->Branch("lep_phiso"  ,   &t_lep_phiso);
  outTree->Branch("lep_chiso"  ,   &t_lep_chiso);
  outTree->Branch("lep_nhiso"  ,   &t_lep_nhiso);
  outTree->Branch("lep_rho"    ,   &t_lep_rho);
  outTree->Branch("lep_pdgId"  ,   &t_lep_pdgId);
  outTree->Branch("lep_charge" ,   &t_lep_charge);
  outTree->Branch("lep_isofull",   &t_lep_isofull);
  outTree->Branch("lep_isofull20", &t_lep_isofull20);
  outTree->Branch("lep_isofull25", &t_lep_isofull25);
  outTree->Branch("lep_isofull30", &t_lep_isofull30);
  outTree->Branch("lep_miniiso",   &t_lep_miniiso);
  outTree->Branch("lepSF",         &t_lepSF);
  outTree->Branch("lepSFUnc",      &t_lepSFUnc);
  outTree->Branch("lepIsoSF",      &t_lepIsoSF);
  outTree->Branch("lepIsoSFUnc",   &t_lepIsoSFUnc);

  //trigger info
  Int_t t_etrig, t_mtrig;
  outTree->Branch("etrig" , &t_etrig , "etrig/I");
  outTree->Branch("mtrig" , &t_mtrig , "mtrig/I");
  outTree->Branch("trigSF" ,   &t_trigSF ,     "trigSF/F");
  outTree->Branch("trigSFUnc" , &t_trigSFUnc , "trigSFUnc/F");


  //b-tag counters
  Int_t t_nbjet_sel, t_nbjet_sel_jecup, t_nbjet_sel_jecdn, t_nbjet_sel_jerup, t_nbjet_sel_jerdn, t_nbjet_sel_bup, t_nbjet_sel_bdn, t_nbjet_sel_udsgup, t_nbjet_sel_udsgdn, t_nbjet_sel_quenchup, t_nbjet_sel_quenchdn;
  outTree->Branch("nbjet_sel"       , &t_nbjet_sel       , "nbjet_sel/I"       );
  outTree->Branch("nbjet_sel_jecup" , &t_nbjet_sel_jecup , "nbjet_sel_jecup/I" );
  outTree->Branch("nbjet_sel_jecdn" , &t_nbjet_sel_jecdn , "nbjet_sel_jecdn/I" );
  outTree->Branch("nbjet_sel_jerup" , &t_nbjet_sel_jerup , "nbjet_sel_jerup/I" );
  outTree->Branch("nbjet_sel_jerdn" , &t_nbjet_sel_jerdn , "nbjet_sel_jerdn/I" );
  outTree->Branch("nbjet_sel_bup"   , &t_nbjet_sel_bup   , "nbjet_sel_bup/I"   );
  outTree->Branch("nbjet_sel_bdn"   , &t_nbjet_sel_bdn   , "nbjet_sel_bdn/I"   );
  outTree->Branch("nbjet_sel_udsgup", &t_nbjet_sel_udsgup, "nbjet_sel_udsgup/I");
  outTree->Branch("nbjet_sel_udsgdn", &t_nbjet_sel_udsgdn, "nbjet_sel_udsgdn/I");
  outTree->Branch("nbjet_sel_quenchup", &t_nbjet_sel_quenchup, "nbjet_sel_quenchup/I");
  outTree->Branch("nbjet_sel_quenchdn", &t_nbjet_sel_quenchdn, "nbjet_sel_quenchdn/I");

  //per-jet kinematics
  std::vector<Float_t> t_jet_pt, t_jet_eta, t_jet_phi, t_jet_mass, t_jet_csvv2,t_jet_deltaQuench;
  std::vector<Float_t> t_jet_matchpt, t_jet_matcheta, t_jet_matchphi, t_jet_matchmass;
  std::vector<Int_t> t_jet_flavor, t_jet_flavorForB;
  outTree->Branch("jet_pt"    , &t_jet_pt    );
  outTree->Branch("jet_eta"   , &t_jet_eta   );
  outTree->Branch("jet_phi"   , &t_jet_phi   );
  outTree->Branch("jet_mass"  , &t_jet_mass  );
  outTree->Branch("jet_deltaQuench"  , &t_jet_deltaQuench );
  outTree->Branch("jet_csvv2" , &t_jet_csvv2 );
  outTree->Branch("jet_genpt"    , &t_jet_matchpt    );
  outTree->Branch("jet_geneta"   , &t_jet_matcheta   );
  outTree->Branch("jet_genphi"   , &t_jet_matchphi   );
  outTree->Branch("jet_genmass"  , &t_jet_matchmass  );
  outTree->Branch("jet_flavor"  , &t_jet_flavor  );
  outTree->Branch("jet_flavorB"  , &t_jet_flavorForB  );


  //prepare to loop over events
  Double_t wgtSum(0);
  std::vector<Double_t> allWgtSum;
  int nEntries = (int)lepTree_p->GetEntries();  
  int entryDiv = ((int)(nEntries/20));    
  cout << inURL << " has " << nEntries << " events to process" << endl;
  if(maxEvents>0) { 
    nEntries=TMath::Min(nEntries,maxEvents); 
    cout << "Number of events to process limited to " << nEntries << endl;
  }

  //first loop : ensure that the ncoll weighting will be normalized
  double ncollWgtNorm(1.0);
  if(isMC && !isPP){
    double ncollSum(0.);
    for(int entry = 0; entry < nEntries; entry++){
      hiTree_p->GetEntry(entry);
      ncollSum+=findNcoll(fForestTree.hiBin);
    }
    if(ncollSum>0) ncollWgtNorm=double(nEntries)/ncollSum;
  }
  
  //loop over events
  for(int entry = 0; entry < nEntries; entry++){
    
    if(entryDiv!=0)if(entry%entryDiv == 0) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    globalTree_p->GetEntry(entry);
    lepTree_p->GetEntry(entry);
    pfCandTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);    
    hltTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    if(rhoTree_p) rhoTree_p->GetEntry(entry);
    if(muHLTObj_p ) muHLTObj_p->GetEntry(entry);
    if(eleHLTObj_p)  eleHLTObj_p->GetEntry(entry);

    //event header
    t_run    = fForestTree.run;
    if(blind && !isMC && t_run>=326887) continue;
    t_lumi   = fForestTree.lumi;
    t_event  = fForestTree.evt;
    t_isData = !isMC;

    //vertex
    t_vx     = fForestTree.vx;
    t_vy     = fForestTree.vy;
    t_vz     = fForestTree.vz;

    //centrality and Ncoll weight
    t_cenbin   = 0.;
    t_ncollWgt = 1.0;
    bool isCentralEvent(false);    
    if(!isPP){
      isCentralEvent=(fForestTree.hiBin<30);
      t_cenbin=0.5*fForestTree.hiBin;
      t_ncollWgt=findNcoll(fForestTree.hiBin)*ncollWgtNorm;
    }
    
    //gen level analysis
    t_weight=1.0;
    t_weight_BRW=1.0;
    t_genFiducial=false;
    t_genlep_id=-1;
    t_genlep_pt=-999.; t_genlep_eta=-999.; t_genlep_phi=-999.; 
    t_genjet_pt.clear();   t_genjet_eta.clear();   t_genjet_phi.clear();   t_genjet_m.clear();
    t_meWeights.clear();
    if(isMC) {

      //gen level selection
      int nlFromTopW(0),nlightlFromTopW(0);
      TLorentzVector topP4(0,0,0,0),antitopP4(0,0,0,0);
      for(size_t i=0; i<fForestGen.mcPID->size(); i++) {
        int pid=fForestGen.mcPID->at(i);
        int mom_pid=fForestGen.mcMomPID->at(i);
        int gmom_pid=fForestGen.mcGMomPID->at(i);
        
        TLorentzVector p4(0,0,0,0);
        p4.SetPtEtaPhiM( fForestGen.mcPt->at(i), fForestGen.mcEta->at(i), fForestGen.mcPhi->at(i), fForestGen.mcMass->at(i) );

        if(pid==6)  topP4=p4;
        if(pid==-6) antitopP4=p4;
        bool isFromTop(abs(mom_pid)==6 || abs(gmom_pid)==6 );
        bool isFromW( abs(mom_pid)==24 || abs(gmom_pid)==24 );
        bool isTauFeedDown( abs(mom_pid)==15 );
            
        //electrons and muons from t->W->l or W->tau->l
        if(abs(pid)==11 || abs(pid)==13) {
          if(isFromTop && isFromW) nlFromTopW++;
          if(isFromW && (isTauFeedDown || isFromTop ) && nlightlFromTopW==0) {
            nlightlFromTopW++;
            if(p4.Pt()>20 && fabs(p4.Eta())<2.4 ) {
              t_genlep_id  = abs(pid);
              t_genlep_pt  = p4.Pt();
              t_genlep_eta = p4.Eta();
              t_genlep_phi = p4.Phi();
            }
          }
        }
        
        //special case for tau leptons (use neutrinos as taus are not stored)
        if(abs(pid)==16 && isFromW && isTauFeedDown) nlFromTopW++;
      }
      
      //gen jet selection
      for (int genjetIter = 0; genjetIter < fForestJets.ngen; genjetIter++){
        TLorentzVector p4;
        p4.SetPtEtaPhiM( fForestJets.genpt[genjetIter],fForestJets.geneta[genjetIter],fForestJets.genphi[genjetIter],fForestJets.genm[genjetIter]);
        if(p4.Pt()<30 || fabs(p4.Eta())>2.0) continue;
        if(nlightlFromTopW) {
          float dR=sqrt( pow(p4.Eta()-t_genlep_eta,2) + pow(TVector2::Phi_mpi_pi(p4.Phi()-t_genlep_phi),2) );
          if(dR<0.4) continue;
        }
        t_genjet_pt.push_back(p4.Pt());
        t_genjet_eta.push_back(p4.Eta());
        t_genjet_phi.push_back(p4.Phi());
        t_genjet_m.push_back(p4.M());
      }

      //consider the event to be in the fiducial volume if it has 1 ch. lepton + >=4 jets
      t_genFiducial=(nlightlFromTopW==1 && (t_genlep_id==11 || t_genlep_id==13) && t_genjet_pt.size()>3); 

      //generator level weights
      t_weight = fForestTree.ttbar_w->size()==0 ? 1. : fForestTree.ttbar_w->at(meIdxList[0]);
      wgtSum += t_weight;    
      if(allWgtSum.size()==0) allWgtSum.resize(meIdxList.size(),0.);
      for(size_t i=0; i<meIdxList.size(); i++) {
        Double_t iwgt(fForestTree.ttbar_w->size()<i  || fForestTree.ttbar_w->size() == 0 ? 1. : fForestTree.ttbar_w->at(meIdxList[i]));
        allWgtSum[i]+=iwgt;
      }

      //BR correction
      t_weight_BRW=isAMCATNLO ? getMadgraphBRWlCorrection(nlFromTopW) : 1.0;

      //top pt modelling (pp-based)
      float topPtWgt = TMath::Exp(0.199-0.00166*topP4.Pt());
      topPtWgt *= TMath::Exp(0.199-0.00166*antitopP4.Pt());
      topPtWgt = TMath::Sqrt(topPtWgt);

      //emulate different top mass from the BW
      std::vector<float> obsm={float(topP4.M()),float(antitopP4.M())};
      float topMassDnWgt = weightBW(rbwigner,obsm,171.5,1.28,172.5,1.31);
      float topMassUpWgt = weightBW(rbwigner,obsm,173.5,1.34,172.5,1.31);
      
      t_meWeights.push_back(topPtWgt);                                                                                                                                
      t_meWeights.push_back(1./topPtWgt);
      t_meWeights.push_back(topMassUpWgt);
      t_meWeights.push_back(topMassDnWgt);
      if(fForestTree.ttbar_w->size()>0) { 
        float nomWgt=fForestTree.ttbar_w->at(meIdxList[0]);
        for(auto meIdx : meIdxList){
          if(meIdx< fForestTree.ttbar_w->size()){
            t_meWeights.push_back( fForestTree.ttbar_w->at(meIdx)/nomWgt );
          }
        }
      }
    }
    
    //start RECO-level event selection
    //global filters
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

    //make a global rho estimate for isolation purposes
    SlimmedPFCollection_t pfColl;
    for(size_t ipf=0; ipf<fForestPF.pfId->size(); ipf++) {
      int id(abs(fForestPF.pfId->at(ipf)));
      float mass(0.13957);  //pions
      if(id==4) mass=0.;    //photons
      if(id>=5) mass=0.497; //K0L
      pfColl.push_back( getSlimmedPF( id, fForestPF.pfPt->at(ipf),fForestPF.pfEta->at(ipf),fForestPF.pfPhi->at(ipf),mass) );
    }
    Float_t globalrho = getRho(pfColl,{1,2,3,4,5,6},-1.,5.);
    t_globalrho = globalrho;

    //lepton selection
    std::vector<LeptonSummary> selMuons,selElectrons;
    
    //muons
    std::vector<TLorentzVector> muHLTP4;
    if(muHLTObjs) muHLTP4=muHLTObjs->getHLTObjectsP4();
    for(unsigned int muIter = 0; muIter < fForestLep.muPt->size(); ++muIter) {
      
      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      float rawpt(fForestLep.muPt->at(muIter));
      p4.SetPtEtaPhiM(rawpt,fForestLep.muEta->at(muIter),fForestLep.muPhi->at(muIter),0.1057);
      if(TMath::Abs(p4.Eta()) > muEtaCut) continue;
      if(p4.Pt() < lepPtCut) continue;

      //require trigger matching
      bool isTrigMatch(false);
      for(auto hp4: muHLTP4) {
        if(hp4.DeltaR(p4)>0.1) continue;
        isTrigMatch=true;
        break;
      }
      if(!isTrigMatch) continue;
      
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

      //save a summary
      LeptonSummary l(13,p4);
      l.rawpt = p4.Pt();
      l.isTrigMatch=isTrigMatch;
      l.charge  = fForestLep.muCharge->at(muIter);
      l.chiso   = fForestLep.muPFChIso->at(muIter);
      l.nhiso   = fForestLep.muPFNeuIso->at(muIter);
      l.phoiso  = fForestLep.muPFPhoIso->at(muIter);
      l.isofull = l.chiso+l.nhiso+l.phoiso;
      int   tmp_rhoind  = getRhoIndex(p4.Eta(),t_etaMin,t_etaMax);
      l.rho = isPP ? globalrho : t_rho->at(tmp_rhoind);
      l.isofullR=getIsolationFull( pfColl, l.p4);
      l.miniiso = getMiniIsolation( pfColl ,l.p4, l.id);
      l.d0      = fForestLep.muD0   ->at(muIter);
      l.d0err   = 0.; //fForestLep.muD0Err->at(muIter); // no d0err for muons!!!
      l.dz      = fForestLep.muDz   ->at(muIter);
      l.origIdx = muIter;
      l.isMatched=false;
      l.isTauFeedDown=false;
      l.idFlags=1;
      selMuons.push_back(l);
    }
   
    //select electrons
    //cf. https://twiki.cern.ch/twiki/pub/CMS/HiHighPt2019/HIN_electrons2018_followUp.pdf
    std::vector<TLorentzVector> eleHLTP4;
    if(eleHLTObjs) eleHLTP4=eleHLTObjs->getHLTObjectsP4() ;
    for(unsigned int eleIter = 0; eleIter < fForestLep.elePt->size(); ++eleIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      float rawpt(fForestLep.elePt->at(eleIter));
      p4.SetPtEtaPhiM(rawpt,fForestLep.eleEta->at(eleIter),fForestLep.elePhi->at(eleIter),0.000511);

      //calibrate pt
      float calpt=calibratedPt(rawpt, p4.Eta(), t_cenbin, isMC);
      p4 *= calpt/rawpt;
      if(TMath::Abs(p4.Eta()) > eleEtaCut) continue;
      if(TMath::Abs(p4.Eta()) > barrelEndcapEta[0] && TMath::Abs(p4.Eta()) < barrelEndcapEta[1] ) continue;
      if(p4.Pt() < lepPtCut) continue;

      //require trigger matching
      bool isTrigMatch(false);
      for(auto hp4: eleHLTP4) {
        if(hp4.DeltaR(p4)>0.2) continue;
        isTrigMatch=true;
        break;
      }
      if(!isTrigMatch) continue;

      //identification
      int idFlags=getElectronId(TMath::Abs(fForestLep.eleSCEta->at(eleIter))< barrelEndcapEta[0],
                                 fForestLep.eleSigmaIEtaIEta->at(eleIter),
                                 fForestLep.eledEtaSeedAtVtx->at(eleIter),
                                 fForestLep.eledPhiAtVtx->at(eleIter),
                                 fForestLep.eleHoverEBc->at(eleIter),
                                 fForestLep.eleEoverPInv->at(eleIter),
                                 fForestLep.eleIP3D->at(eleIter),
                                 fForestLep.eleMissHits->at(eleIter),
                                 isCentralEvent);
      
      //id'ed electron
      if(!isLooseElectron(idFlags)) continue;

      //save a summary
      LeptonSummary l(11,p4);
      l.rawpt = p4.Pt()*(rawpt/calpt);
      l.idFlags=idFlags;
      l.isTrigMatch=isTrigMatch;
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
      l.isofull = l.chiso+l.nhiso+l.phoiso;
      int   tmp_rhoind  = getRhoIndex(p4.Eta(),t_etaMin,t_etaMax);
      l.rho = isPP ? globalrho : t_rho->at(tmp_rhoind);
      l.isofullR= getIsolationFull( pfColl, l.p4);
      l.miniiso = getMiniIsolation( pfColl ,l.p4, l.id);
      l.d0      = fForestLep.eleD0   ->at(eleIter);
      l.d0err   = fForestLep.eleD0Err->at(eleIter);
      l.dz      = fForestLep.eleDz   ->at(eleIter);
      l.origIdx=eleIter;
      l.isMatched=false;
      l.isTauFeedDown=false;
      selElectrons.push_back(l);
    }
    
    //sort selected leptons by pt and require at least one
    std::sort(selMuons.begin(),     selMuons.end(),     orderByPt);
    std::sort(selElectrons.begin(), selElectrons.end(), orderByPt);
    if(selMuons.size()+selElectrons.size()==0) continue;

    //set the final lepton candidate giving prefence to muons, then electrons
    //in case of ambiguity use electron only if pT is higher
    int selCh(-1);
    LeptonSummary *selLepton=0;
    if(selMuons.size()>0 && mtrig) {
      selCh=13;
      selLepton=&(selMuons[0]);
    }
    if(selElectrons.size()>0 && etrig) {
      if(selMuons.size()==0 || (selCh==13 && selElectrons[0].p4.Pt()>selMuons[0].p4.Pt())) {
        selCh=11;
        selLepton=&(selElectrons[0]);
      }
    }
    t_etrig  = etrig;
    t_mtrig  = mtrig;
    if(selLepton==0) continue;

    //lepton info
    t_lep_pt    = selLepton->p4.Pt();
    t_lep_eta   = selLepton->p4.Eta() ;
    t_lep_phi   = selLepton->p4.Phi() ;
    t_lep_idflags=selLepton->idFlags;
    t_lep_d0    = selLepton->d0 ;
    t_lep_d0err = selLepton->d0err;
    t_lep_dz    = selLepton->dz  ;
    t_lep_chiso = selLepton->chiso ;
    t_lep_phiso = selLepton->phoiso ;
    t_lep_nhiso = selLepton->nhiso ;
    t_lep_rho   = selLepton->rho ;
    t_lep_pdgId = selLepton->id ;
    t_lep_charge= selLepton->charge ;
    t_lep_isofull= selLepton->isofull ;
    t_lep_isofull20= selLepton->isofullR[0] ;
    t_lep_isofull25= selLepton->isofullR[1] ;
    t_lep_isofull30= selLepton->isofullR[2] ;
    t_lep_miniiso = selLepton->miniiso ;

    //get expected trigger efficiencies and measured scale factors
    std::pair<float,float> ltrigEff, ltrigSF;
    float pt(selLepton->p4.Pt()),eta(selLepton->p4.Eta()),abseta(fabs(eta));
    if(selCh==11){
      ltrigEff=std::pair<float,float>(e_mctrigeff->Eval(pt),0.0);
      ltrigSF=eleEff.eval(pt, abseta<barrelEndcapEta[0], t_cenbin, true, false); //HLT (L1 is unity by definition in this trigger menu)
    }else{
      ltrigEff=std::pair<float,float>(tnp_weight_trig_pbpb(pt,eta,300),0.0);
      float mutrigSF=tnp_weight_trig_pbpb(pt,eta,0);
      float deltaTnp=max(fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,-1)),fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,-2)));
      deltaTnp += pow(0.05,2);                                                             //HTL centrality dependence
      float deltaStat=max(fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,1)),fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,2)));
      float deltaSF=sqrt(deltaTnp*deltaTnp+deltaStat*deltaStat);
      ltrigSF=std::pair<float,float>(mutrigSF,deltaSF);        
    }
    
    t_trigSF     = ltrigSF.first;
    t_trigSFUnc  = ltrigSF.second;

    //reco/tracking+id scale factors
    float sfVal(1.0),sfValUnc(0.0);
    if(selCh==13) {
      sfVal=tnp_weight_muid_pbpb( selLepton->p4.Eta(), 0 );                         //central value
      sfValUnc += pow(fabs(tnp_weight_muid_pbpb( selLepton->p4.Eta(),+1)-sfVal),2); //stat
      sfValUnc += pow(fabs(tnp_weight_muid_pbpb( selLepton->p4.Eta(),-1)-sfVal),2); //syst
      sfValUnc += pow(0.006,2);                                                           //tracking uncertainty
      sfValUnc += pow(0.02,2);                                                            //tracking centrality dependence
      sfValUnc += pow(0.01,2);                                                            //identification centrality dependence
      sfValUnc = sqrt(sfValUnc);
    }else {
      std::pair<float,float > eleIDsf=eleEff.eval(selLepton->p4.Pt(), fabs(selLepton->p4.Eta())<barrelEndcapEta[0], t_cenbin, false, false); //ID
      sfVal=eleIDsf.first;
      sfValUnc=eleIDsf.second;
      std::pair<float,float > eleRECOsf=eleEff.eval(selLepton->p4.Pt(), fabs(selLepton->p4.Eta())<barrelEndcapEta[0], t_cenbin, false, true); //RECO 
      sfValUnc = sqrt(pow(sfValUnc/sfVal,2)+pow(eleRECOsf.second/eleRECOsf.first,2));
      sfVal*=eleRECOsf.first;
    }    
    t_lepSF=sfVal;
    t_lepSFUnc=sfValUnc;

    //isolation efficiencies
    TString isoKey(isCentralEvent ? "cen" : "periph");
    isoKey+=abs(selLepton->id)==13 ? "_169" : "_121";
    Int_t xbin=isoEffSFs[isoKey]->GetXaxis()->FindBin( min(selLepton->p4.Pt(), isoEffSFs[isoKey]->GetXaxis()->GetXmax()-0.01) );
    Int_t ybin=isoEffSFs[isoKey]->GetYaxis()->FindBin( min(fabs(selLepton->p4.Eta()), isoEffSFs[isoKey]->GetYaxis()->GetXmax()) );
    if (isoEffSFs[isoKey]->GetBinContent(xbin,ybin) != 0.){
      t_lepIsoSF   = isoEffSFs[isoKey]->GetBinContent(xbin,ybin) ;
      t_lepIsoSFUnc= isoEffSFs[isoKey]->GetBinError  (xbin,ybin) ;
    } else{
      t_lepIsoSF   = 1.00;
      t_lepIsoSFUnc= 0.05;
    }

    //skip bad runs, no trigger run
    if(isSingleMuPD || isMuSkimedMCPD) {
      if(std::find(badMuonTriggerRuns.begin(), badMuonTriggerRuns.end(), fForestTree.run) != badMuonTriggerRuns.end() and !isMuSkimedMCPD) continue;
      if(mtrig==0) continue;
    }
    if(isSingleElePD || isEleSkimedMCPD) {
      if(etrig==0) continue;
    }

    //analyze jets
    int npfjets(0); 
    std::vector<BtagInfo_t> pfJetsIdx;
    std::vector<TLorentzVector> pfJetsP4;
    t_nbjet_sel = 0; 
    t_nbjet_sel_jecup    = 0; t_nbjet_sel_jecdn    = 0;
    t_nbjet_sel_jerup    = 0; t_nbjet_sel_jerdn    = 0;
    t_nbjet_sel_bup      = 0; t_nbjet_sel_bdn      = 0;
    t_nbjet_sel_udsgup   = 0; t_nbjet_sel_udsgdn   = 0;
    t_nbjet_sel_quenchup = 0; t_nbjet_sel_quenchdn = 0;
    t_jet_pt   .clear();
    t_jet_eta  .clear();
    t_jet_phi  .clear();
    t_jet_mass .clear();
    t_jet_deltaQuench.clear();
    t_jet_csvv2.clear();
    t_jet_matchpt  .clear();
    t_jet_matcheta .clear();
    t_jet_matchphi .clear();
    t_jet_matchmass.clear();
    t_jet_flavor.clear();
    t_jet_flavorForB.clear();
    for(int jetIter = 0; jetIter < fForestJets.nref; jetIter++){

      //at least two tracks
      if(fForestJets.trackN[jetIter]<2) continue;

      TLorentzVector jp4(0,0,0,0);
      if(isMC){
	JECMC.SetJetPT(fForestJets.rawpt[jetIter]);
	JECMC.SetJetEta(fForestJets.jteta[jetIter]);
	JECMC.SetJetPhi(fForestJets.jtphi[jetIter]);
	jp4.SetPtEtaPhiM( JECMC.GetCorrectedPT(),fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);
	JEUMC.SetJetPT(JECMC.GetCorrectedPT());
	JEUMC.SetJetEta(fForestJets.jteta[jetIter]);
	JEUMC.SetJetPhi(fForestJets.jtphi[jetIter]);
      }
      else{
	JECData.SetJetPT(fForestJets.rawpt[jetIter]);
	JECData.SetJetEta(fForestJets.jteta[jetIter]);
	JECData.SetJetPhi(fForestJets.jtphi[jetIter]);
	jp4.SetPtEtaPhiM( JECData.GetCorrectedPT(),fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);
	JEUData.SetJetPT(JECMC.GetCorrectedPT());
	JEUData.SetJetEta(fForestJets.jteta[jetIter]);
	JEUData.SetJetPhi(fForestJets.jtphi[jetIter]);
      }

      float csvVal=fForestJets.discr_csvV2[jetIter];
      int nsvtxTk=fForestJets.svtxntrk[jetIter];
      float msvtx=fForestJets.svtxm[jetIter];

      if(jp4.Pt()<20.) continue; // smaller pT cut here to avoid the full loop
      if(fabs(jp4.Eta())>2.0) continue;

      //cross clean wrt to selected lepton
      if(jp4.DeltaR(selLepton->p4)<0.4) continue;

      bool isBTagged(csvVal>csvWPList[csvWP]);      
      
      // simple matching to the closest jet in dR. require at least dR < 0.3
      TLorentzVector matchjp4(0,0,0,0);
      int refFlavor(0),refFlavorForB(0);
      if (isMC){
        for (int genjetIter = 0; genjetIter < fForestJets.ngen; genjetIter++){
          if (jetIter == fForestJets.genmatchindex[genjetIter]) {
            matchjp4.SetPtEtaPhiM( fForestJets.genpt[genjetIter],fForestJets.geneta[genjetIter],fForestJets.genphi[genjetIter],fForestJets.genm[genjetIter]);
          }
        }  
        refFlavor=fForestJets.refparton_flavor[jetIter];
        refFlavorForB=fForestJets.refparton_flavorForB[jetIter];
      }
      
      pfJetsIdx.push_back(std::make_tuple(pfJetsP4.size(),nsvtxTk,msvtx,csvVal,matchjp4,refFlavor,refFlavorForB));
      pfJetsP4.push_back(jp4);
      npfjets++;
            
      //b-tag counting variations
      if (jp4.Pt() > 30. && isBTagged) t_nbjet_sel       += 1;
      if (isMC){
        if (jp4.Pt() > 30 * (1 + JEUMC.GetUncertainty().first) && isBTagged)  t_nbjet_sel_jecup += 1;
        if (jp4.Pt() > 30 * (1 - JEUMC.GetUncertainty().second) && isBTagged) t_nbjet_sel_jecdn += 1;
        float cjer(0.);
        if ( abs(refFlavorForB) ) cjer = 1. + (1.2 -1.) * (jp4.Pt() - matchjp4.Pt()) / jp4.Pt(); // hard coded 1.2
        else cjer = smearRand->Gaus(1., 0.2);
        if (jp4.Pt()*cjer > 30. && isBTagged) t_nbjet_sel_jerup += 1;
        if (jp4.Pt()/cjer > 30. && isBTagged) t_nbjet_sel_jerdn += 1;
        quenchingModel->SetParameter(0, 50.); // this sets the omega_c parameter. if we want to make this centrality dependent
        float tmp_quench_loss = quenchingModel->GetRandom();
        tmp_quench_loss = TMath::Abs(TMath::Sin(jp4.Theta())*tmp_quench_loss); // make it only on the transverse part...
        // make it centrality dependent
        float centralitySuppression = centralityModel->Eval(t_cenbin/100.);
        if (jp4.Pt()                  > 30. && isBTagged)                       t_nbjet_sel_quenchup += 1;
        if (jp4.Pt()-tmp_quench_loss*centralitySuppression  > 30. && isBTagged) t_nbjet_sel_quenchdn += 1;

        bool isBTaggedNew(0);
        float tmp_btageff = btagEfficiencies(refFlavorForB, t_cenbin);
        if (jp4.Pt() < 30.) continue;
        t_nbjet_sel_bup    = t_nbjet_sel; t_nbjet_sel_bdn    = t_nbjet_sel;
        t_nbjet_sel_udsgup = t_nbjet_sel; t_nbjet_sel_udsgdn = t_nbjet_sel;

        if ( abs(refFlavorForB) == 5){

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 1.05, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_bup = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_bup = t_nbjet_sel-1;

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 0.95, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_bdn = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_bdn = t_nbjet_sel-1;
        }
        else {

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 1.15, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_udsgup = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_udsgup = t_nbjet_sel-1;

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 0.85, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_udsgdn = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_udsgdn = t_nbjet_sel-1;
        }
      }
    }
    // set the bjet variations to the nominal one for data
    if (!isMC){
      t_nbjet_sel_jecup    = t_nbjet_sel; t_nbjet_sel_jecdn    = t_nbjet_sel;
      t_nbjet_sel_jerup    = t_nbjet_sel; t_nbjet_sel_jerdn    = t_nbjet_sel;
      t_nbjet_sel_bup      = t_nbjet_sel; t_nbjet_sel_bdn      = t_nbjet_sel;
      t_nbjet_sel_udsgup   = t_nbjet_sel; t_nbjet_sel_udsgdn   = t_nbjet_sel;
      t_nbjet_sel_quenchup = t_nbjet_sel; t_nbjet_sel_quenchdn = t_nbjet_sel;
    }
    if(npfjets<4) continue;

    //order jets by their b-tag info content and save them
    std::sort(pfJetsIdx.begin(),       pfJetsIdx.end(),      orderByBtagInfo);
    for (size_t ij = 0; ij < pfJetsIdx.size(); ij++) {
      int idx = std::get<0>(pfJetsIdx[ij]);
      t_jet_pt   .push_back( pfJetsP4[idx].Pt()  );
      t_jet_eta  .push_back( pfJetsP4[idx].Eta() );
      t_jet_phi  .push_back( pfJetsP4[idx].Phi() );
      t_jet_mass .push_back( pfJetsP4[idx].M()   );
      t_jet_csvv2.push_back( std::get<3>(pfJetsIdx[ij])   );      
      t_jet_matchpt  .push_back( std::get<4>(pfJetsIdx[ij]).Pt());
      t_jet_matcheta .push_back( std::get<4>(pfJetsIdx[ij]).Eta());
      t_jet_matchphi .push_back( std::get<4>(pfJetsIdx[ij]).Phi());
      t_jet_matchmass.push_back( std::get<4>(pfJetsIdx[ij]).M());
      t_jet_flavor.push_back( std::get<5>(pfJetsIdx[ij]) );
      t_jet_flavorForB.push_back( std::get<6>(pfJetsIdx[ij]) );


      quenchingModel->SetParameter(0, 50.); // this sets the omega_c parameter. if we want to make this centrality dependent
      float tmp_quench_loss = quenchingModel->GetRandom();
      tmp_quench_loss = TMath::Abs(TMath::Sin(pfJetsP4[idx].Theta())*tmp_quench_loss); // make it only on the transverse part...
      // make it centrality dependent
      float centralitySuppression = centralityModel->Eval(t_cenbin/100.);
      t_jet_deltaQuench.push_back( tmp_quench_loss*centralitySuppression);
    }
        
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

    fOut->Close();
  }

  return 0;
}
