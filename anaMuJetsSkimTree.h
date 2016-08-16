#ifndef anaMuJetsSkimTree_h
#define anaMuJetsSkimTree_h

#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>



class anaMuJetsSkimTree{

public:

  TTree* tr_[5];     //Change, or think of something, when more files are included

  anaMuJetsSkimTree();
  ~anaMuJetsSkimTree(); 

  //Interface functions
  void SetCSVCut(double CSVcut) { CSVCut_ = CSVcut; };
  void SetnumbB(double numbB) { numbB_ = numbB; };
  
  double GetCSVCut() {return CSVCut_; };
  int GethiBin() {return hiBin; };
  double GetjetPt(int index) {return jetPt[index]; };
  double GetlepPt(int index) {return lepPt[index]; };
  double GetlepIso(int index) {return lepIso[index]; };
  int GetnLep() {return nLep; };

  void OpenFiles(TString FileName, int i);
  void BuildHistograms();
  int FindLeadingMuonAfterCuts(int option);
  std::vector<int> FindJetsAfterCuts(int indexMuon);
  std::vector<int> FindbJets(std::vector<int> indexJets);
  void CalculateNormalizationHistograms(TH1F* h[4], int option);
  void NormalizeHistograms(int option);
  void NormalizeHistogramsToOne(int i, int nBjets);
  double CalculateDeltaPhi(int ilep, int index_b);
  void SetAddressBranches(int i);
  std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets);
  std::vector<int> BuildpTVectorLeadingbJets(std::vector<int> indexbJets);
  void FillHistogramsBeforeCuts(int option);
  void FillHistogramsAfterCuts(int option);

private:
  //Constant variables
  const Float_t muM = .1056583715;

  //Objects
  TFile* f_[5];      //Change, or think of something, when more files are included

  //Global variables
  unsigned int run_, lumi_;
  ULong64_t evt_;
  int hiBin;
  float vz;
  
  //Max leptons in skimFile = 2, to be sure I used 4
  int nLep;
  int lepID[4];
  float lepPt[4];
  float lepPhi[4];
  float lepEta[4];
  int lepChg[4];
  float lepIso[4];

  //Max jets = 500
  int nJt;
  double sum_jtPt;
  int numb_bJets;
  float jetPt[500];
  float jtPhi[500];
  float jtEta[500];
  float jtM[500];
  float discr_csvV1[500];


  //Cut values
  int numbB_ = 1;
  double CSVCut_ = 0.75;
  double mupTCut_ = 18;
  double jtpTCut_ = 30.;
  double jtEtaCut_ = 2.;
  double drJetToMuonCut_ = 0.3;
};

#endif