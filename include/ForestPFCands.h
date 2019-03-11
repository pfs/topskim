#ifndef ForestPFCands_h
#define ForestPFCands_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

class ForestPFCands {
public :
  ForestPFCands(TChain *t)
    {
      t->SetBranchAddress("pfId", &pfId);
      t->SetBranchAddress("pfPt", &pfPt);
      t->SetBranchAddress("pfEta", &pfEta);
      t->SetBranchAddress("pfPhi", &pfPhi);
      t->SetBranchAddress("pfM", &pfM);
      t->SetBranchAddress("trkAlgo", &trkAlgo);
      t->SetBranchAddress("trkPtError", &trkPtError);
      t->SetBranchAddress("trkNHit", &trkNHit);
      t->SetBranchAddress("trkChi2", &trkChi2);
      t->SetBranchAddress("trkNdof", &trkNdof);
    }
  ~ForestPFCands() {}

   std::vector<int>     *pfId =0;
   std::vector<float>   *pfPt =0;
   std::vector<float>   *pfEta =0;
   std::vector<float>   *pfPhi =0;
   std::vector<float>   *pfM =0;
   std::vector<int>     *trkAlgo =0;
   std::vector<float>   *trkPtError =0;
   std::vector<float>   *trkNHit =0;
   std::vector<float>   *trkChi2 =0;
   std::vector<float>   *trkNdof =0;
};



#endif 
