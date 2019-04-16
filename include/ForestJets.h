//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  8 23:52:51 2019 by ROOT version 6.12/07
// from TTree t/akPu4PFpatJetsWithBtagging Jet Analysis Tree
// found on file: /eos/cms/store/group/phys_top/PbPbTTbar_2018/SkimElectrons/Chunk_0_ext0.root
//////////////////////////////////////////////////////////

#ifndef ForestJets_h
#define ForestJets_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLeaf.h>

// Header file for the classes stored in the TTree if any.

class ForestJets {
public :
  ForestJets(TChain *t) {
    t->SetBranchAddress("nref", &nref);
    t->SetBranchAddress("rawpt", rawpt);
    t->SetBranchAddress("jtpt", jtpt);
    t->SetBranchAddress("jteta", jteta);
    t->SetBranchAddress("jty", jty);
    t->SetBranchAddress("jtphi", jtphi);
    t->SetBranchAddress("jtpu", jtpu);
    t->SetBranchAddress("jtm", jtm);
    t->SetBranchAddress("jtarea", jtarea);
    t->SetBranchAddress("jtPfCHF", jtPfCHF);
    t->SetBranchAddress("jtPfNHF", jtPfNHF);
    t->SetBranchAddress("jtPfCEF", jtPfCEF);
    t->SetBranchAddress("jtPfNEF", jtPfNEF);
    t->SetBranchAddress("jtPfMUF", jtPfMUF);
    t->SetBranchAddress("jtPfCHM", jtPfCHM);
    t->SetBranchAddress("jtPfNHM", jtPfNHM);
    t->SetBranchAddress("jtPfCEM", jtPfCEM);
    t->SetBranchAddress("jtPfNEM", jtPfNEM);
    t->SetBranchAddress("jtPfMUM", jtPfMUM);
    t->SetBranchAddress("jttau1", jttau1);
    t->SetBranchAddress("jttau2", jttau2);
    t->SetBranchAddress("jttau3", jttau3);
    t->SetBranchAddress("discr_jetID_cuts", discr_jetID_cuts);
    t->SetBranchAddress("discr_jetID_bdt", discr_jetID_bdt);
    t->SetBranchAddress("discr_fr01", discr_fr01);
    t->SetBranchAddress("trackMax", trackMax);
    t->SetBranchAddress("trackSum", trackSum);
    t->SetBranchAddress("trackN", trackN);
    t->SetBranchAddress("trackHardSum", trackHardSum);
    t->SetBranchAddress("trackHardN", trackHardN);
    t->SetBranchAddress("chargedMax", chargedMax);
    t->SetBranchAddress("chargedSum", chargedSum);
    t->SetBranchAddress("chargedN", chargedN);
    t->SetBranchAddress("chargedHardSum", chargedHardSum);
    t->SetBranchAddress("chargedHardN", chargedHardN);
    t->SetBranchAddress("photonMax", photonMax);
    t->SetBranchAddress("photonSum", photonSum);
    t->SetBranchAddress("photonN", photonN);
    t->SetBranchAddress("photonHardSum", photonHardSum);
    t->SetBranchAddress("photonHardN", photonHardN);
    t->SetBranchAddress("neutralMax", neutralMax);
    t->SetBranchAddress("neutralSum", neutralSum);
    t->SetBranchAddress("neutralN", neutralN);
    t->SetBranchAddress("hcalSum", hcalSum);
    t->SetBranchAddress("ecalSum", ecalSum);
    t->SetBranchAddress("eMax", eMax);
    t->SetBranchAddress("eSum", eSum);
    t->SetBranchAddress("eN", eN);
    t->SetBranchAddress("muMax", muMax);
    t->SetBranchAddress("muSum", muSum);
    t->SetBranchAddress("muN", muN);
    t->SetBranchAddress("discr_ssvHighEff", discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur", discr_ssvHighPur);
    t->SetBranchAddress("discr_csvV1", discr_csvV1);
    t->SetBranchAddress("discr_csvV2", discr_csvV2);
    t->SetBranchAddress("discr_muByIp3", discr_muByIp3);
    t->SetBranchAddress("discr_muByPt", discr_muByPt);
    t->SetBranchAddress("discr_prob", discr_prob);
    t->SetBranchAddress("discr_probb", discr_probb);
    t->SetBranchAddress("discr_tcHighEff", discr_tcHighEff);
    t->SetBranchAddress("discr_tcHighPur", discr_tcHighPur);
    t->SetBranchAddress("ndiscr_ssvHighEff", ndiscr_ssvHighEff);
    t->SetBranchAddress("ndiscr_ssvHighPur", ndiscr_ssvHighPur);
    t->SetBranchAddress("ndiscr_csvV1", ndiscr_csvV1);
    t->SetBranchAddress("ndiscr_csvV2", ndiscr_csvV2);
    t->SetBranchAddress("ndiscr_muByPt", ndiscr_muByPt);
    t->SetBranchAddress("pdiscr_csvV1", pdiscr_csvV1);
    t->SetBranchAddress("pdiscr_csvV2", pdiscr_csvV2);
    t->SetBranchAddress("nsvtx", nsvtx);
    t->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
    t->SetBranchAddress("refparton_flavor", refparton_flavor);
    TLeaf *leaf = t->GetLeaf("svtxntrk");
    if(leaf) {
      TString lname(leaf->GetTypeName());
      if(lname=="Int_t") {
        t->SetBranchAddress("svtxntrk", svtxntrk);
        t->SetBranchAddress("svtxdl", svtxdl);
        t->SetBranchAddress("svtxdls", svtxdls);
        t->SetBranchAddress("svtxdl2d", svtxdl2d);
        t->SetBranchAddress("svtxdls2d", svtxdls2d);
        t->SetBranchAddress("svtxm", svtxm);
        t->SetBranchAddress("svtxpt", svtxpt);
      }else {
        /*
          t->SetBranchAddress("svType", &svtype_vec);
          t->SetBranchAddress("svtxntrk", &svtxntrk_vec);
          t->SetBranchAddress("svtxdl", &svtxdl_vec);
          t->SetBranchAddress("svtxdls", &svtxdls_vec);
          t->SetBranchAddress("svtxdl2d", &svtxdl2d_vec);
          t->SetBranchAddress("svtxdls2d", &svtxdls2d_vec);
          t->SetBranchAddress("svtxm", &svtxm_vec);
          t->SetBranchAddress("svtxpt", &svtxpt_vec);
        */
      }
    }

    t->SetBranchAddress("svtxmcorr", svtxmcorr);
    t->SetBranchAddress("nIPtrk", nIPtrk);
    t->SetBranchAddress("nselIPtrk", nselIPtrk);
    t->SetBranchAddress("mue", mue);
    t->SetBranchAddress("mupt", mupt);
    t->SetBranchAddress("mueta", mueta);
    t->SetBranchAddress("muphi", muphi);
    t->SetBranchAddress("mudr", mudr);
    t->SetBranchAddress("muptrel", muptrel);
    t->SetBranchAddress("muchg", muchg);
    t->SetBranchAddress("ngen"          , &ngen         );
    t->SetBranchAddress("genpt"         , genpt        );
    t->SetBranchAddress("geneta"        , geneta       );
    t->SetBranchAddress("genphi"        , genphi       );
    t->SetBranchAddress("genm"          , genm         );
    t->SetBranchAddress("genmatchindex" , genmatchindex);

  }
  ~ForestJets() {}


   // Declaration of leaf types
   Int_t           nref;
   Float_t         rawpt[85];   //[nref]
   Float_t         jtpt[85];   //[nref]
   Float_t         jteta[85];   //[nref]
   Float_t         jty[85];   //[nref]
   Float_t         jtphi[85];   //[nref]
   Float_t         jtpu[85];   //[nref]
   Float_t         jtm[85];   //[nref]
   Float_t         jtarea[85];   //[nref]
   Float_t         jtPfCHF[85];   //[nref]
   Float_t         jtPfNHF[85];   //[nref]
   Float_t         jtPfCEF[85];   //[nref]
   Float_t         jtPfNEF[85];   //[nref]
   Float_t         jtPfMUF[85];   //[nref]
   Int_t           jtPfCHM[85];   //[nref]
   Int_t           jtPfNHM[85];   //[nref]
   Int_t           jtPfCEM[85];   //[nref]
   Int_t           jtPfNEM[85];   //[nref]
   Int_t           jtPfMUM[85];   //[nref]
   Float_t         jttau1[85];   //[nref]
   Float_t         jttau2[85];   //[nref]
   Float_t         jttau3[85];   //[nref]
   Float_t         discr_jetID_cuts[85];   //[nref]
   Float_t         discr_jetID_bdt[85];   //[nref]
   Float_t         discr_fr01[85];   //[nref]
   Float_t         trackMax[85];   //[nref]
   Float_t         trackSum[85];   //[nref]
   Int_t           trackN[85];   //[nref]
   Float_t         trackHardSum[85];   //[nref]
   Int_t           trackHardN[85];   //[nref]
   Float_t         chargedMax[85];   //[nref]
   Float_t         chargedSum[85];   //[nref]
   Int_t           chargedN[85];   //[nref]
   Float_t         chargedHardSum[85];   //[nref]
   Int_t           chargedHardN[85];   //[nref]
   Float_t         photonMax[85];   //[nref]
   Float_t         photonSum[85];   //[nref]
   Int_t           photonN[85];   //[nref]
   Float_t         photonHardSum[85];   //[nref]
   Int_t           photonHardN[85];   //[nref]
   Float_t         neutralMax[85];   //[nref]
   Float_t         neutralSum[85];   //[nref]
   Int_t           neutralN[85];   //[nref]
   Float_t         hcalSum[85];   //[nref]
   Float_t         ecalSum[85];   //[nref]
   Float_t         eMax[85];   //[nref]
   Float_t         eSum[85];   //[nref]
   Int_t           eN[85];   //[nref]
   Float_t         muMax[85];   //[nref]
   Float_t         muSum[85];   //[nref]
   Int_t           muN[85];   //[nref]
   Float_t         discr_ssvHighEff[85];   //[nref]
   Float_t         discr_ssvHighPur[85];   //[nref]
   Float_t         discr_csvV1[85];   //[nref]
   Float_t         discr_csvV2[85];   //[nref]
   Float_t         discr_muByIp3[85];   //[nref]
   Float_t         discr_muByPt[85];   //[nref]
   Float_t         discr_prob[85];   //[nref]
   Float_t         discr_probb[85];   //[nref]
   Float_t         discr_tcHighEff[85];   //[nref]
   Float_t         discr_tcHighPur[85];   //[nref]
   Float_t         ndiscr_ssvHighEff[85];   //[nref]
   Float_t         ndiscr_ssvHighPur[85];   //[nref]
   Float_t         ndiscr_csvV1[85];   //[nref]
   Float_t         ndiscr_csvV2[85];   //[nref]
   Float_t         ndiscr_muByPt[85];   //[nref]
   Float_t         pdiscr_csvV1[85];   //[nref]
   Float_t         pdiscr_csvV2[85];   //[nref]
   Int_t           nsvtx[85];   //[nref]
   Int_t           svtxntrk[85];   //[nref]
   Float_t         svtxdl[85];   //[nref]
   Float_t         svtxdls[85];   //[nref]
   Float_t         svtxdl2d[85];   //[nref]
   Float_t         svtxdls2d[85];   //[nref]
   Float_t         svtxm[85];   //[nref]
   Float_t         svtxpt[85];   //[nref]
   Int_t refparton_flavorForB[85],refparton_flavor[85];
   std::vector<std::vector<int> > *svtype_vec;
   std::vector<std::vector<int> > *svtxntrk_vec;
   std::vector<std::vector<float> > *svtxdl_vec;
   std::vector<std::vector<float> > *svtxdls_vec;
   std::vector<std::vector<float> > *svtxdl2d_vec;
   std::vector<std::vector<float> > *svtxdls2d_vec;
   std::vector<std::vector<float> > *svtxm_vec;
   std::vector<std::vector<float> > *svtxpt_vec;
   Float_t         svtxmcorr[85];   //[nref]
   Int_t           nIPtrk[85];   //[nref]
   Int_t           nselIPtrk[85];   //[nref]
   Float_t         mue[85];   //[nref]
   Float_t         mupt[85];   //[nref]
   Float_t         mueta[85];   //[nref]
   Float_t         muphi[85];   //[nref]
   Float_t         mudr[85];   //[nref]
   Float_t         muptrel[85];   //[nref]
   Int_t           muchg[85];   //[nref]
   Int_t           ngen;
   Float_t         genpt[85];
   Float_t         geneta[85];
   Float_t         genphi[85];
   Float_t         genm[85];
   Int_t           genmatchindex[85];
};

#endif
