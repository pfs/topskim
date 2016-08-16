#ifndef histControlDistributions_h
#define histControlDistributions_h

#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

  TH1F* h_events_[5];
  TH1F* h_totevents[5];

  TH1F* h_lepPt_[5];
  TH1F* h_lepPseu_[5];
  TH1F* h_jetHt_[5];
  TH1F* h_jetCSV_[5];
  TH1F* h_lepbjetMinv_min_[5];
  TH1F* h_ttMinv_min_[5];
  TH1F* h_qqMinv_min_[5];
  TH1F* h_pTselectbjet_[5];
  TH1F* h_Phi_1stb_[5];
  TH1F* h_Phi_allb_[5];
  TH1F* h_Phi_b_minpi_[5];
  TH2F* h_2dCSV_[5];
  TH2F* h_2dCSV2_[5];
  TGraph* CSV_effcut;
  TGraph* CSV_effcut2;

  THStack* lepPt_Stack;
  THStack* lepPseu_Stack;
  THStack* jetHt_Stack;
  THStack* jetCSV_Stack;
  THStack* lepbjetMinv_min_Stack;
  THStack* ttMinv_min_Stack;
  THStack* pTselectbjet_Stack;
  THStack* qqMinv_min_Stack;

#endif