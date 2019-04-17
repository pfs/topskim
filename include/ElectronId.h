#ifndef _electronid_h_
#define _electronid_h_

#include "TMath.h"
#include <vector>

//electron id (separate for EB and EE, depending on centrality)
//cf. https://twiki.cern.ch/twiki/pub/CMS/HiEgamma2019/HIN_electronID_17042019_EGamma.pdf
unsigned int getElectronId(bool isEB,float sihih, float detaVtx, float dphiVtx, float hoe, float eop, float d0, float dz, float missHits,bool isCentralEvent) {

  std::vector<float> max_sihih(4,0.),max_detaVtx(4,0.),max_dphiVtx(4,0.), max_hoe(4,0.), max_eop(4,0.), max_d0(4,0.), max_dz(4,0.), max_misshits(4,0.);
  if(isEB){
    if(isCentralEvent){
      max_sihih    = {0.0164, 0.0161, 0.0161, 0.0155};
      max_detaVtx  = {0.0065, 0.0053, 0.0053, 0.0044};
      max_dphiVtx  = {0.0610, 0.0288, 0.0147, 0.0094};
      max_hoe      = {0.3470, 0.1984, 0.1735, 0.1569};
      max_eop      = {0.1157, 0.1129, 0.0098, 0.0047};
      max_misshits = {3,      1,      1,      1};
      max_d0       = {999.,   999.,   999.,   999. };
      max_dz       = {999.,   999.,   999.,   999. };
    }
    else{
      max_sihih    = {0.0140, 0.0117, 0.0111, 0.0110};
      max_detaVtx  = {0.0081, 0.0071, 0.0041, 0.0026};
      max_dphiVtx  = {0.0302, 0.0221, 0.0157, 0.0107};
      max_hoe      = {0.1993, 0.1892, 0.1428, 0.1367};
      max_eop      = {0.0589, 0.0405, 0.0064, 0.0063};
      max_misshits = {3,      1,      1,      1};
      max_d0       = {999.,   999.,   999.,   999. };
      max_dz       = {999.,   999.,   999.,   999. };
    }
  }
  else{
    if(isCentralEvent){
      max_sihih    = {0.0481, 0.0479, 0.0477, 0.0450};
      max_detaVtx  = {0.0146, 0.0145, 0.0090, 0.0090};
      max_dphiVtx  = {0.1025, 0.0516, 0.0334, 0.0174};
      max_hoe      = {0.2838, 0.1910, 0.1824, 0.1820};
      max_eop      = {0.0129, 0.0115, 0.0047, 0.0040};
      max_misshits = {3,      1,      1,      1};
      max_d0       = {0.02,   0.02,   0.02,   0.02};
      max_dz       = {0.04,   0.04,   0.04,   0.04};
    }
    else{
      max_sihih    = {0.0460, 0.0447, 0.0440, 0.0422};
      max_detaVtx  = {0.0108, 0.0108, 0.0075, 0.0066};
      max_dphiVtx  = {0.0576, 0.0301, 0.0191, 0.0126};
      max_hoe      = {0.2104, 0.1627, 0.1616, 0.1473};
      max_eop      = {0.0754, 0.0281, 0.0086, 0.0037};
      max_misshits = {3,      1,      1,      1};
      max_d0       = {0.02,   0.02,   0.02,   0.02};
      max_dz       = {0.04,   0.04,   0.04,   0.04};
    }
  }

  //set the bits according to the ids passed
  unsigned int idFlags(0);
  for(size_t i=0; i<4; i++) {
    bool passId(sihih<max_sihih[i] &&
                fabs(detaVtx)<max_detaVtx[i] &&
                fabs(dphiVtx)<max_dphiVtx[i] &&
                hoe<max_hoe[i] &&
                eop<max_eop[i] &&
                missHits<=max_misshits[i] &&
                fabs(d0)<max_d0[i] &&
                fabs(dz)<max_dz[i] );
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
