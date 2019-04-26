#ifndef ForestSkim_h
#define ForestSkim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

class ForestSkim {
public :
 ForestSkim(TChain *t)
      {      
	t->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
	t->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
    t->SetBranchAddress("phfCoincFilter2Th4", &phfCoincFilter2Th4);
    t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
	t->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	t->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);

      }

  ~ForestSkim() {}
  // Declaration of leaf types                                                                                                                                                                             
  int phfCoincFilter;
  int HBHENoiseFilterResult;
  int pprimaryVertexFilter;
  int pcollisionEventSelection;
  int phfCoincFilter2Th4;
  int pclusterCompatibilityFilter;
};



#endif 
