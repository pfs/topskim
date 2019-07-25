#ifndef _electronid_h_
#define _electronid_h_

#include "TMath.h"
#include <vector>

//electron id (separate for EB and EE, depending on centrality)
//cf. https://twiki.cern.ch/twiki/pub/CMS/HiEgamma2019/HIN_electronID_17042019_EGamma.pdf
unsigned int getElectronId(bool isEB,float sihih, float detaseedVtx, float dphiVtx, float hoebc, float eop, float ip3d,  float missHits, bool isCentralEvent) {

  std::vector<float> max_sihih(4,0.),max_detaseedVtx(4,0.),max_dphiVtx(4,0.), max_hoebc(4,0.), max_eop(4,0.), max_ip3d(4,0.),  max_misshits(4,0.);
  if(isEB){
    if(isCentralEvent){
      max_sihih    = {0.0147, 0.0135, 0.0116, 0.0104};
      max_detaseedVtx  = {0.0041, 0.0038, 0.0037, 0.0029};
      max_dphiVtx  = {0.0853, 0.0376, 0.0224, 0.0206};
      max_hoebc      = {0.2733, 0.1616, 0.1589, 0.1459};
      max_eop      = {0.0367, 0.0177, 0.0173, 0.0105};
      max_misshits = {3,      1,      1,      1};
      max_ip3d       = {0.03,   0.03,   0.03,   0.03};
    }
    else{
      max_sihih    = {0.0113, 0.0107, 0.0101, 0.0099};
      max_detaseedVtx  = {0.0037, 0.0035, 0.0033, 0.0026};
      max_dphiVtx  = {0.1280, 0.0327, 0.0210, 0.0170};
      max_hoebc      = {0.1814, 0.1268, 0.0311, 0.0067};
      max_eop      = {0.1065, 0.0774, 0.0701, 0.0077};
      max_misshits = {3,      1,      1,      1};
      max_ip3d       = {0.03,   0.03,   0.03,   0.03};
    }
  }
  else{
    if(isCentralEvent){
      max_sihih    = {0.048, 0.0466, 0.0418, 0.0358};
      max_detaseedVtx  = {0.0097, 0.0063, 0.0062, 0.0051};
      max_dphiVtx  = {0.2348, 0.1186, 0.0373, 0.0266};
      max_hoebc      = {0.1898, 0.1317, 0.1092, 0.0925};
      max_eop      = {0.0300, 0.0201, 0.0133, 0.0065};
      max_misshits = {3,      1,      1,      1};
      max_ip3d       = {0.03,   0.03,   0.03,   0.03};
    }
    else{
      max_sihih    = {0.0376, 0.0339, 0.0316, 0.0288};
      max_detaseedVtx  = {0.0074, 0.0067, 0.0051, 0.0044};
      max_dphiVtx  = {0.2085, 0.0838, 0.0384, 0.0266};
      max_hoebc      = {0.1138, 0.0977, 0.0810, 0.0655};
      max_eop      = {0.0237, 0.0193, 0.0192, 0.0123};
      max_misshits = {3,      1,      1,      1};
      max_ip3d       = {0.03,   0.03,   0.03,   0.03};
    }
  }

  //set the bits according to the ids passed
  unsigned int idFlags(0);
  for(size_t i=0; i<4; i++) {
    bool passId(sihih<max_sihih[i] &&
                fabs(detaseedVtx)<max_detaseedVtx[i] &&
                fabs(dphiVtx)<max_dphiVtx[i] &&
                hoebc<max_hoebc[i] &&
                eop<max_eop[i] &&
                missHits<=max_misshits[i] &&
                fabs(ip3d)<max_ip3d[i] 
                );
    idFlags |= (passId<<i);
  }

  return idFlags;
}

//
bool isVetoElectron(unsigned int idFlags)   { return (idFlags & 0x1);      }
bool isLooseElectron(unsigned int idFlags)  { return ((idFlags>>1) & 0x1); }
bool isMediumElectron(unsigned int idFlags) { return ((idFlags>>2) & 0x1); }
bool isTightElectron(unsigned int idFlags)  { return ((idFlags>>3) & 0x1); }

#endif
