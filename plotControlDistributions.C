#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"
#include "histControlDistributions.h"

#include <string>
#include <vector>

void LayoutHistograms();
TCanvas* PlotHistograms(int nbJets);

void plotControlDistributions(std::string outFileName = "", int nJets = 4, int nbJets = 1){

  anaMuJetsSkimTree ana;

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/Data.root",0);			//LepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_tt/MCtt.root", 1);			//NoLepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_W/MCW.root", 2);			//NoLepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_DY/MCDY.root", 3);			//NoLepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/DataMultijets.root", 4);	//NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  //cout << "Making Control Distribution plots in the following pT interval." << endl;
  //cout << "     ptLow = " << pTLow << " , pTHigh = " << pTHigh << endl;
  cout << endl << "Making Control Distribution plots for " << nJets << "j" << nbJets << "b." << endl;

  if(nbJets == 0){

  	cout << "   So we are not placing a cut on the CSVv1 variable" << endl << endl;
  	ana.BuildHistograms();
  	for(int i = 0; i < 5; i++)  ana.FillHistogramsBeforeCuts(i);

  } else {
  	
  	cout << "   We are placing a cut on the CSVv1 variable to tag b-jets" << endl << endl;
  	ana.BuildHistograms();
  	for(int i = 0; i < 5; i++)  ana.FillHistogramsAfterCuts(i);

  }

  LayoutHistograms();
  for(int i = 0; i < 5; i++)  ana.NormalizeHistograms(i);
  //for(int i = 0; i < 4; i++) NormalizeHistogramsToOne(i, nbjets);
  TCanvas* cst = PlotHistograms(nbJets);

  cout << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  outFile->Write();
  cst->Write();
  delete outFile;
}




void LayoutHistograms(){

  Color_t colors[4] = {920,800,42,32};

  h_lepPt_[0]->SetMarkerStyle(20);
  h_lepPseu_[0]->SetMarkerStyle(20);
  h_jetHt_[0]->SetMarkerStyle(20);
  h_jetCSV_[0]->SetMarkerStyle(20);
  h_lepbjetMinv_min_[0]->SetMarkerStyle(20);
  h_ttMinv_min_[0]->SetMarkerStyle(20);
  h_qqMinv_min_[0]->SetMarkerStyle(20);
  h_pTselectbjet_[0]->SetMarkerStyle(20);

  for(int i = 1; i < 5; i++){
  	h_lepPt_[i]->SetFillColor(colors[i-1]);
    h_lepPseu_[i]->SetFillColor(colors[i-1]);
    h_jetHt_[i]->SetFillColor(colors[i-1]);
    h_jetCSV_[i]->SetFillColor(colors[i-1]);
    h_lepbjetMinv_min_[i]->SetFillColor(colors[i-1]);
    h_ttMinv_min_[i]->SetFillColor(colors[i-1]);
    h_qqMinv_min_[i]->SetFillColor(colors[i-1]);
    h_pTselectbjet_[i]->SetFillColor(colors[i-1]);
  }
}


TCanvas* PlotHistograms(int nbJets){

  TCanvas *cst = new TCanvas("cst","Histograms",10,10,1800,1000);
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  lepPt_Stack = new THStack("lepPt_Stack","pT distributions (muon) for 1l+X (at least 4j1b);pT;Counts");
  lepPseu_Stack = new THStack("lepPseu_Stack","Pseudorapidity distributions (muon) for 1l+X (at least 4j1b);|eta|;Counts");
  jetHt_Stack = new THStack("jetHt_Stack","H_T distributions for 1l+X (at least 4j1b);H_T;Counts");
  jetCSV_Stack = new THStack("jetCSV_Stack","CSV distributions for all jets in 1l+X (at least 4j0b);CSVv1;Counts");
  pTselectbjet_Stack = new THStack("pTselectbjet_Stack","pT distributions for tagged b-jet in 1l+X (at least 4j1b);pT;Counts");
  lepbjetMinv_min_Stack = new THStack("lepbjetMinv_min_Stack","M_lb distributions (at least 4j1b);M_lb;Counts");
  ttMinv_min_Stack = new THStack("ttMinv_min_Stack","M_tt distributions (at least 4j1b);M_tt;Counts");
  qqMinv_min_Stack = new THStack("qqMinv_min_Stack","M_qq distributions (at least 4j1b);M_qq;Counts");

  for(int j = 3; j > 0; j--){
    lepPt_Stack->Add(h_lepPt_[j]);
    lepPseu_Stack->Add(h_lepPseu_[j]);
    jetHt_Stack->Add(h_jetHt_[j]);
    
    if(nbJets == 0) jetCSV_Stack->Add(h_jetCSV_[j]);
    else { 
      pTselectbjet_Stack->Add(h_pTselectbjet_[j]);
      lepbjetMinv_min_Stack->Add(h_lepbjetMinv_min_[j]);
      ttMinv_min_Stack->Add(h_ttMinv_min_[j]);
      qqMinv_min_Stack->Add(h_qqMinv_min_[j]);
    }
  }
  lepPt_Stack->Add(h_lepPt_[4]);
  lepPseu_Stack->Add(h_lepPseu_[4]);
  jetHt_Stack->Add(h_jetHt_[4]);
  if(nbJets == 0) jetCSV_Stack->Add(h_jetCSV_[4]);
  else { 
  	pTselectbjet_Stack->Add(h_pTselectbjet_[4]);
    lepbjetMinv_min_Stack->Add(h_lepbjetMinv_min_[4]);
    ttMinv_min_Stack->Add(h_ttMinv_min_[4]);
    qqMinv_min_Stack->Add(h_qqMinv_min_[4]);
  }

  leg->AddEntry(h_lepPt_[0],"Data","p");
  leg->AddEntry(h_lepPt_[1],"ttbar (MC)","f");
  leg->AddEntry(h_lepPt_[2],"Wjets (MC)","f");
  leg->AddEntry(h_lepPt_[3],"DY (MC)","f");
  leg->AddEntry(h_lepPt_[4],"Multijets (Data)","f");

  cst->Divide(3,2);

  cst->cd(1);
  lepPt_Stack->Draw(); 
  if(lepPt_Stack->GetMaximum() < h_lepPt_[0]->GetMaximum()){
    lepPt_Stack->SetMaximum((int)( h_lepPt_[0]->GetMaximum() + 0.15 * h_lepPt_[0]->GetMaximum() ));
  }
  h_lepPt_[0]->Draw("same,ep");
  leg->Draw();
    
  cst->cd(2);
  lepPseu_Stack->Draw();
  if(lepPseu_Stack->GetMaximum() < h_lepPseu_[0]->GetMaximum()){
    lepPseu_Stack->SetMaximum((int)( h_lepPseu_[0]->GetMaximum() + 0.15 * h_lepPseu_[0]->GetMaximum() ));
  }
  h_lepPseu_[0]->Draw("same,ep");
  leg->Draw();
  
  cst->cd(3);
  jetHt_Stack->Draw();
  if(jetHt_Stack->GetMaximum() < h_jetHt_[0]->GetMaximum()){
    jetHt_Stack->SetMaximum((int)( h_jetHt_[0]->GetMaximum() + 0.15 * h_jetHt_[0]->GetMaximum() ));
  }
  h_jetHt_[0]->Draw("same,ep");
  leg->Draw();
  
  cst->cd(4);
  if(nbJets == 0){
    jetCSV_Stack->Draw();
    if(jetCSV_Stack->GetMaximum() < h_jetCSV_[0]->GetMaximum()){
      jetCSV_Stack->SetMaximum((int)( h_jetCSV_[0]->GetMaximum() + 0.15 * h_jetCSV_[0]->GetMaximum() ));
    }
    h_jetCSV_[0]->Draw("same,ep");
    leg->Draw();
  } 
  else {
    lepbjetMinv_min_Stack->Draw();
    if(lepbjetMinv_min_Stack->GetMaximum() < h_lepbjetMinv_min_[0]->GetMaximum()){
      lepbjetMinv_min_Stack->SetMaximum((int)( h_lepbjetMinv_min_[0]->GetMaximum() + 0.15 * h_lepbjetMinv_min_[0]->GetMaximum() ));
    }
    h_lepbjetMinv_min_[0]->Draw("same,ep");
    leg->Draw();
    
    cst->cd(5);
    pTselectbjet_Stack->Draw();
    if(pTselectbjet_Stack->GetMaximum() < h_pTselectbjet_[0]->GetMaximum()){
      pTselectbjet_Stack->SetMaximum((int)( h_pTselectbjet_[0]->GetMaximum() + 0.15 * h_pTselectbjet_[0]->GetMaximum() ));
    }
    h_pTselectbjet_[0]->Draw("same,ep");
    leg->Draw();

    cst->cd(6);
    ttMinv_min_Stack->Draw();
    if(ttMinv_min_Stack->GetMaximum() < h_ttMinv_min_[0]->GetMaximum()){
      ttMinv_min_Stack->SetMaximum((int)( h_ttMinv_min_[0]->GetMaximum() + 0.15 * h_ttMinv_min_[0]->GetMaximum() ));
    }
    h_ttMinv_min_[0]->Draw("same,ep");
    leg->Draw();

    cst->cd(7);
    qqMinv_min_Stack->Draw();
    if(qqMinv_min_Stack->GetMaximum() < h_qqMinv_min_[0]->GetMaximum()){
      qqMinv_min_Stack->SetMaximum((int)( h_qqMinv_min_[0]->GetMaximum() + 0.15 * h_qqMinv_min_[0]->GetMaximum() ));
    }
    h_qqMinv_min_[0]->Draw("same,ep");
    leg->Draw();
  }
    
  cst->Update();

  return cst;
}