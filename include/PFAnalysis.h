#ifndef PFAnalysis_h
#define PFAnalysis_h

#include <algorithm>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "TLorentzVector.h"

typedef std::pair<int,TLorentzVector> SlimmedPF_t;
typedef std::vector<SlimmedPF_t> SlimmedPFCollection_t;

//wrapper to return a slimmed PF
SlimmedPF_t getSlimmedPF(int id, float pt, float eta, float phi, float m) {
  TLorentzVector p4;
  p4.SetPtEtaPhiM(pt,eta,phi,m);
  return SlimmedPF_t(id,p4); 
}

//compute FastJet rho for a set of particles in a given pt/eta range
float getRho(SlimmedPFCollection_t &coll, std::vector<int> ids,float minAbsEta=-1,float maxAbsEta=2.4,float minPt=0.5) {

  std::vector<fastjet::PseudoJet> cands;
  for(auto spf : coll) {
    if(spf.second.Pt()<minPt) continue;
    float abseta(spf.second.Eta());
    if(abseta<minAbsEta) continue;
    if(abseta>maxAbsEta) continue;
    if( std::find(ids.begin(),ids.end(),spf.first)==ids.end() ) continue;
    fastjet::PseudoJet ip(spf.second.Px(),spf.second.Py(),spf.second.Pz(),spf.second.E());
    ip.set_user_index(spf.first);
  }

  fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm,0.5);
  fastjet::Selector sel_rap( fastjet::SelectorAbsRapRange(minAbsEta,maxAbsEta) );
  fastjet::AreaDefinition area_def(fastjet::active_area,fastjet::GhostedAreaSpec(minAbsEta,maxAbsEta));
  fastjet::JetMedianBackgroundEstimator jmbe(sel_rap, jet_def_for_rho, area_def);
  jmbe.set_particles(cands);
  
  return jmbe.rho();
}

#endif
