#ifndef ForestHLTObject_h
#define ForestHLTObject_h

#include <TROOT.h>
#include <TChain.h>
#include <vector>

/**
   @short a wrapper to read HLT Object kinematics event-by-event
 */

class ForestHLTObject {
 public :
 ForestHLTObject(TChain *t) :
    isvalid(false){
    
    if(t){
      t->SetBranchAddress("TriggerObjID", &TriggerObjId);
      t->SetBranchAddress("pt", &pt);
      t->SetBranchAddress("eta", &eta);
      t->SetBranchAddress("phi", &phi);
      t->SetBranchAddress("mass", &mass);
      isvalid=true;
    }
  }
  
  std::vector<TLorentzVector> getHLTObjectsP4(){
    std::vector<TLorentzVector> toRet;
    if(!isvalid && pt!=NULL) return toRet;
    for(size_t i=0; i<pt->size(); i++){
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM( pt->at(i), eta->at(i), phi->at(i), mass->at(i) );
      toRet.push_back(p4);
    }
    return toRet;
  }

  bool isValid() { return isvalid; }

  ~ForestHLTObject() { }


 private:
  std::vector<Double_t> *TriggerObjId=0, *pt=0, *eta=0, *phi=0, *mass=0;
  bool isvalid;
};

#endif
