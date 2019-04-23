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
    float abseta( fabs(spf.second.Eta()) );
    if(abseta<minAbsEta) continue;
    if(abseta>maxAbsEta) continue;
    if(ids.size() && std::find(ids.begin(),ids.end(),spf.first)==ids.end() ) continue;
    fastjet::PseudoJet ip(spf.second.Px(),spf.second.Py(),spf.second.Pz(),spf.second.E());
    ip.set_user_index(spf.first);
    cands.push_back(ip);
  }

  fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm,0.5);
  fastjet::Selector sel_rap( fastjet::SelectorAbsRapRange(minAbsEta,maxAbsEta) );
  fastjet::AreaDefinition area_def(fastjet::active_area,fastjet::GhostedAreaSpec(minAbsEta,maxAbsEta));
  fastjet::JetMedianBackgroundEstimator jmbe(sel_rap, jet_def_for_rho, area_def);
  jmbe.set_particles(cands);
  
  return jmbe.rho();
}


//
float getMiniIsolation(SlimmedPFCollection_t &pfCands,
                       TLorentzVector p4, int lid,
                       float r_iso_min=0.05, float r_iso_max=0.2, float kt_scale=6.,
                       bool charged_only=false) 
{

  if (p4.Pt()<5.) return 99999.;
  
  float deadcone_nh(0.015), deadcone_ch(0.015), deadcone_ph(0.015);
  
  //loop over PF candidates to build isolation components
  float iso_nh(0.), iso_ch(0.), iso_ph(0.);
  float ptThresh(0.5);
  if(abs(lid)==11) ptThresh=0;
  float r_iso = (float)TMath::Max((float)r_iso_min,
                                  (float)TMath::Min((float)r_iso_max, (float)(kt_scale/p4.Pt())));  
  for(auto pfc : pfCands) {           

    int pfid(abs(pfc.first));
    float pfpt(pfc.second.Pt());
    if(pfpt<ptThresh) continue;

    float dr = p4.DeltaR(pfc.second);
    if (dr > r_iso) continue;
    
    //photons
    if(pfid==4) {
      if( dr < deadcone_ph) continue;
      iso_ph += pfpt;
    }
    
    //neutral hadrons
    if(pfid==5 || pfid==6) {
      if( dr < deadcone_nh) continue;
      iso_nh += pfpt;
    }
    
    if(pfid==1 || pfid==2 || pfid==3) {
      if( dr < deadcone_ch) continue;
      iso_ch += pfpt;
    }
  }

  //sum up isolation components
  float iso(0.);
  if (charged_only){
    iso = iso_ch;
  } else {
    iso = iso_ph+iso_nh+iso_ch;      
  }

  return iso/p4.Pt();
}


//
std::vector<float> getIsolationFull(SlimmedPFCollection_t &pfCands,
                                    TLorentzVector p4, 
                                    std::vector<float> rMax={0.2,0.25,0.3},
                                    float deadCone=0.015,
                                    float ptThresh=0.5)
{

  std::vector<float> iso(rMax.size(),0.);
  
  if (p4.Pt()<5.) return iso;
  for(auto pfc : pfCands) {           

    float pfpt(pfc.second.Pt());
    if(pfpt<ptThresh) continue;
    
    for(size_t i=0; i<rMax.size(); i++) {
      float dr = p4.DeltaR(pfc.second);
      if (dr > rMax[i]) continue;
      if (dr< deadCone) continue;
      iso[i] += pfpt;
    }
  } 
  
  return iso;
}


#endif
