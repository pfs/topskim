#ifndef ForestGen_h
#define ForestGen_h

#include <iostream>
#include <vector>

class ForestGen {
public :
 ForestGen(TChain *t) 
   : 
  nMC(0),
    mcPID(0),
    mcStatus(0),
    mcMomPID(0),
    mcPt(0),
    mcEta(0),
    mcPhi(0),
    mcMass(0)
     {

       t->SetBranchStatus("nMC", 1);
       t->SetBranchStatus("mc*", 1);
       t->SetBranchAddress("nMC", &nMC);
       t->SetBranchAddress("mcPID", &mcPID);
       t->SetBranchAddress("mcStatus", &mcStatus);
       t->SetBranchAddress("mcMomPID", &mcMomPID);
       t->SetBranchAddress("mcPt", &mcPt);
       t->SetBranchAddress("mcEta", &mcEta);       
       t->SetBranchAddress("mcPhi", &mcPhi);
       t->SetBranchAddress("mcMass", &mcMass);
     };
   ~ForestGen(){};

   //gen info
   Int_t           nMC;
   std::vector<Int_t>     *mcPID,*mcStatus,*mcMomPID;
   std::vector<Float_t>   *mcPt,*mcEta,*mcPhi,*mcMass;
};
#endif
