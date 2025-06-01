//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 16 18:22:34 2018 by ROOT version 6.06/01
// from TTree hits/HGC rechits
// found on file: muon_v10.root
//////////////////////////////////////////////////////////

#ifndef HGCNtupleVariables_h
#define HGCNtupleVariables_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class HGCNtupleVariables {
public:
  HGCNtupleVariables(TTree * /*tree*/ = 0) : fChain(0) {}
  ~HGCNtupleVariables() {}
  // void    Init(TTree *tree);
  void Init(TTree *tree, TTree *tree2);
  Bool_t Notify();
  Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }

  TTree *fChain;   //! pointer to the analyzed TTree or TChain
  TTree *fChain2;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //! current Tree number in a TChain
  Int_t fCurrent2; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  //
   vector<float>   *iEtaEle1;
   vector<float>   *iPhiEle1;
   vector<float>   *Hit_ES_Eta_Ele1;
   vector<float>   *Hit_ES_Phi_Ele1;
   vector<float>   *Hit_ES_X_Ele1;
   vector<float>   *Hit_ES_Y_Ele1;
   vector<float>   *Hit_ES_Z_Ele1;
   vector<float>   *ES_RecHitEnEle1;
   vector<float>   *Hit_Eta_Ele1;
   vector<float>   *Hit_Phi_Ele1;
   vector<float>   *Hit_X_Ele1;
   vector<float>   *Hit_Y_Ele1;
   vector<float>   *Hit_Z_Ele1;
   vector<float>   *RecHitEnEle1;
   vector<float>   *RecHitFracPho1;
   vector<int>     *RecHitGain1;
   vector<bool>    *RecHitQuality1;
   vector<float>   *HitNoiseEle1;
   vector<float>   *iEtaEle2;
   vector<float>   *iPhiEle2;
   vector<float>   *Hit_ES_Eta_Ele2;
   vector<float>   *Hit_ES_Phi_Ele2;
   vector<float>   *Hit_ES_X_Ele2;
   vector<float>   *Hit_ES_Y_Ele2;
   vector<float>   *Hit_ES_Z_Ele2;
   vector<float>   *ES_RecHitEnEle2;
   vector<float>   *Hit_Eta_Ele2;
   vector<float>   *Hit_Phi_Ele2;
   vector<float>   *Hit_X_Ele2;
   vector<float>   *Hit_Y_Ele2;
   vector<float>   *Hit_Z_Ele2;
   vector<float>   *RecHitEnEle2;
   vector<float>   *RecHitFracPho2;
   vector<int>     *RecHitGain2;
   vector<bool>    *RecHitQuality2;
   vector<float>   *HitNoiseEle2;
   vector<float>   *iEtaEle3;
   vector<float>   *iPhiEle3;
   vector<float>   *Hit_ES_Eta_Ele3;
   vector<float>   *Hit_ES_Phi_Ele3;
   vector<float>   *Hit_ES_X_Ele3;
   vector<float>   *Hit_ES_Y_Ele3;
   vector<float>   *Hit_ES_Z_Ele3;
   vector<float>   *ES_RecHitEnEle3;
   vector<float>   *Hit_Eta_Ele3;
   vector<float>   *Hit_Phi_Ele3;
   vector<float>   *Hit_X_Ele3;
   vector<float>   *Hit_Y_Ele3;
   vector<float>   *Hit_Z_Ele3;
   vector<float>   *RecHitEnEle3;
   vector<float>   *RecHitFracPho3;
   vector<int>     *RecHitGain3;
   vector<bool>    *RecHitQuality3;
   vector<float>   *HitNoiseEle3;
   vector<float>   *iEtaEle4;
   vector<float>   *iPhiEle4;
   vector<float>   *Hit_ES_Eta_Ele4;
   vector<float>   *Hit_ES_Phi_Ele4;
   vector<float>   *Hit_ES_X_Ele4;
   vector<float>   *Hit_ES_Y_Ele4;
   vector<float>   *Hit_ES_Z_Ele4;
   vector<float>   *ES_RecHitEnEle4;
   vector<float>   *Hit_Eta_Ele4;
   vector<float>   *Hit_Phi_Ele4;
   vector<float>   *Hit_X_Ele4;
   vector<float>   *Hit_Y_Ele4;
   vector<float>   *Hit_Z_Ele4;
   vector<float>   *RecHitEnEle4;
   vector<float>   *RecHitFracPho4;
   vector<int>     *RecHitGain4;
   vector<bool>    *RecHitQuality4;
   vector<float>   *HitNoiseEle4;
   Int_t           nElectrons;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *energy;
   vector<float>   *energy_ecal_mustache;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
   vector<int>     *passMVAMediumId;
   vector<float>   *Ele_R9;
   vector<float>   *Ele_S4;
   vector<float>   *Ele_SigIEIE;
   vector<float>   *Ele_SigIPhiIPhi;
   vector<float>   *Ele_SCEtaW;
   vector<float>   *Ele_SCPhiW;
   vector<float>   *Ele_CovIEtaIEta;
   vector<float>   *Ele_CovIEtaIPhi;
   vector<float>   *Ele_ESSigRR;
   vector<float>   *Ele_SCRawE;
   vector<float>   *Ele_SC_ESEnByRawE;
   vector<float>   *Ele_HadOverEm;
   vector<float>   *Ele_Gen_Pt;
   vector<float>   *Ele_Gen_Eta;
   vector<float>   *Ele_Gen_Phi;
   vector<float>   *Ele_Gen_E;
   Float_t         rho;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;

// List of branches
   TBranch        *b_iEtaEle1;   //!
   TBranch        *b_iPhiEle1;   //!
   TBranch        *b_Hit_ES_Eta_Ele1;   //!
   TBranch        *b_Hit_ES_Phi_Ele1;   //!
   TBranch        *b_Hit_ES_X_Ele1;   //!
   TBranch        *b_Hit_ES_Y_Ele1;   //!
   TBranch        *b_Hit_ES_Z_Ele1;   //!
   TBranch        *b_ES_RecHitEnEle1;   //!
   TBranch        *b_Hit_Eta_Ele1;   //!
   TBranch        *b_Hit_Phi_Ele1;   //!
   TBranch        *b_Hit_X_Ele1;   //!
   TBranch        *b_Hit_Y_Ele1;   //!
   TBranch        *b_Hit_Z_Ele1;   //!
   TBranch        *b_RecHitEnEle1;   //!
   TBranch        *b_RecHitFracPho1;   //!
   TBranch        *b_RecHitGain1;   //!
   TBranch        *b_RecHitQuality1;   //!
   TBranch        *b_HitNoiseEle1;   //!
   TBranch        *b_iEtaEle2;   //!
   TBranch        *b_iPhiEle2;   //!
   TBranch        *b_Hit_ES_Eta_Ele2;   //!
   TBranch        *b_Hit_ES_Phi_Ele2;   //!
   TBranch        *b_Hit_ES_X_Ele2;   //!
   TBranch        *b_Hit_ES_Y_Ele2;   //!
   TBranch        *b_Hit_ES_Z_Ele2;   //!
   TBranch        *b_ES_RecHitEnEle2;   //!
   TBranch        *b_Hit_Eta_Ele2;   //!
   TBranch        *b_Hit_Phi_Ele2;   //!
   TBranch        *b_Hit_X_Ele2;   //!
   TBranch        *b_Hit_Y_Ele2;   //!
   TBranch        *b_Hit_Z_Ele2;   //!
   TBranch        *b_RecHitEnEle2;   //!
   TBranch        *b_RecHitFracPho2;   //!
   TBranch        *b_RecHitGain2;   //!
   TBranch        *b_RecHitQuality2;   //!
   TBranch        *b_HitNoiseEle2;   //!
   TBranch        *b_iEtaEle3;   //!
   TBranch        *b_iPhiEle3;   //!
   TBranch        *b_Hit_ES_Eta_Ele3;   //!
   TBranch        *b_Hit_ES_Phi_Ele3;   //!
   TBranch        *b_Hit_ES_X_Ele3;   //!
   TBranch        *b_Hit_ES_Y_Ele3;   //!
   TBranch        *b_Hit_ES_Z_Ele3;   //!
   TBranch        *b_ES_RecHitEnEle3;   //!
   TBranch        *b_Hit_Eta_Ele3;   //!
   TBranch        *b_Hit_Phi_Ele3;   //!
   TBranch        *b_Hit_X_Ele3;   //!
   TBranch        *b_Hit_Y_Ele3;   //!
   TBranch        *b_Hit_Z_Ele3;   //!
   TBranch        *b_RecHitEnEle3;   //!
   TBranch        *b_RecHitFracPho3;   //!
   TBranch        *b_RecHitGain3;   //!
   TBranch        *b_RecHitQuality3;   //!
   TBranch        *b_HitNoiseEle3;   //!
   TBranch        *b_iEtaEle4;   //!
   TBranch        *b_iPhiEle4;   //!
   TBranch        *b_Hit_ES_Eta_Ele4;   //!
   TBranch        *b_Hit_ES_Phi_Ele4;   //!
   TBranch        *b_Hit_ES_X_Ele4;   //!
   TBranch        *b_Hit_ES_Y_Ele4;   //!
   TBranch        *b_Hit_ES_Z_Ele4;   //!
   TBranch        *b_ES_RecHitEnEle4;   //!
   TBranch        *b_Hit_Eta_Ele4;   //!
   TBranch        *b_Hit_Phi_Ele4;   //!
   TBranch        *b_Hit_X_Ele4;   //!
   TBranch        *b_Hit_Y_Ele4;   //!
   TBranch        *b_Hit_Z_Ele4;   //!
   TBranch        *b_RecHitEnEle4;   //!
   TBranch        *b_RecHitFracPho4;   //!
   TBranch        *b_RecHitGain4;   //!
   TBranch        *b_RecHitQuality4;   //!
   TBranch        *b_HitNoiseEle4;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_energy_ecal_mustache;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_passTightId;   //!
   TBranch        *b_passMVAMediumId;   //!
   TBranch        *b_Ele_R9;   //!
   TBranch        *b_Ele_S4;   //!
   TBranch        *b_Ele_SigIEIE;   //!
   TBranch        *b_Ele_SigIPhiIPhi;   //!
   TBranch        *b_Ele_SCEtaW;   //!
   TBranch        *b_Ele_SCPhiW;   //!
   TBranch        *b_Ele_CovIEtaIEta;   //!
   TBranch        *b_Ele_CovIEtaIPhi;   //!
   TBranch        *b_Ele_ESSigRR;   //!
   TBranch        *b_Ele_SCRawE;   //!
   TBranch        *b_Ele_SC_ESEnByRawE;   //!
   TBranch        *b_Ele_HadOverEm;   //!
   TBranch        *b_Ele_Gen_Pt;   //!
   TBranch        *b_Ele_Gen_Eta;   //!
   TBranch        *b_Ele_Gen_Phi;   //!
   TBranch        *b_Ele_Gen_E;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
}; //Modified by Somanko

#endif

#ifdef HGCNtupleVariables_cxx

void HGCNtupleVariables::Init(TTree *tree, TTree *tree2) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
   iEtaEle1 = 0;
   iPhiEle1 = 0;
   Hit_ES_Eta_Ele1 = 0;
   Hit_ES_Phi_Ele1 = 0;
   Hit_ES_X_Ele1 = 0;
   Hit_ES_Y_Ele1 = 0;
   Hit_ES_Z_Ele1 = 0;
   ES_RecHitEnEle1 = 0;
   Hit_Eta_Ele1 = 0;
   Hit_Phi_Ele1 = 0;
   Hit_X_Ele1 = 0;
   Hit_Y_Ele1 = 0;
   Hit_Z_Ele1 = 0;
   RecHitEnEle1 = 0;
   RecHitFracPho1 = 0;
   RecHitGain1 = 0;
   RecHitQuality1 = 0;
   HitNoiseEle1 = 0;
   iEtaEle2 = 0;
   iPhiEle2 = 0;
   Hit_ES_Eta_Ele2 = 0;
   Hit_ES_Phi_Ele2 = 0;
   Hit_ES_X_Ele2 = 0;
   Hit_ES_Y_Ele2 = 0;
   Hit_ES_Z_Ele2 = 0;
   ES_RecHitEnEle2 = 0;
   Hit_Eta_Ele2 = 0;
   Hit_Phi_Ele2 = 0;
   Hit_X_Ele2 = 0;
   Hit_Y_Ele2 = 0;
   Hit_Z_Ele2 = 0;
   RecHitEnEle2 = 0;
   RecHitFracPho2 = 0;
   RecHitGain2 = 0;
   RecHitQuality2 = 0;
   HitNoiseEle2 = 0;
   iEtaEle3 = 0;
   iPhiEle3 = 0;
   Hit_ES_Eta_Ele3 = 0;
   Hit_ES_Phi_Ele3 = 0;
   Hit_ES_X_Ele3 = 0;
   Hit_ES_Y_Ele3 = 0;
   Hit_ES_Z_Ele3 = 0;
   ES_RecHitEnEle3 = 0;
   Hit_Eta_Ele3 = 0;
   Hit_Phi_Ele3 = 0;
   Hit_X_Ele3 = 0;
   Hit_Y_Ele3 = 0;
   Hit_Z_Ele3 = 0;
   RecHitEnEle3 = 0;
   RecHitFracPho3 = 0;
   RecHitGain3 = 0;
   RecHitQuality3 = 0;
   HitNoiseEle3 = 0;
   iEtaEle4 = 0;
   iPhiEle4 = 0;
   Hit_ES_Eta_Ele4 = 0;
   Hit_ES_Phi_Ele4 = 0;
   Hit_ES_X_Ele4 = 0;
   Hit_ES_Y_Ele4 = 0;
   Hit_ES_Z_Ele4 = 0;
   ES_RecHitEnEle4 = 0;
   Hit_Eta_Ele4 = 0;
   Hit_Phi_Ele4 = 0;
   Hit_X_Ele4 = 0;
   Hit_Y_Ele4 = 0;
   Hit_Z_Ele4 = 0;
   RecHitEnEle4 = 0;
   RecHitFracPho4 = 0;
   RecHitGain4 = 0;
   RecHitQuality4 = 0;
   HitNoiseEle4 = 0;
   pt = 0;
   eta = 0;
   phi = 0;
   energy = 0;
   energy_ecal_mustache = 0;
   passMediumId = 0;
   passTightId = 0;
   passMVAMediumId = 0;
   Ele_R9 = 0;
   Ele_S4 = 0;
   Ele_SigIEIE = 0;
   Ele_SigIPhiIPhi = 0;
   Ele_SCEtaW = 0;
   Ele_SCPhiW = 0;
   Ele_CovIEtaIEta = 0;
   Ele_CovIEtaIPhi = 0;
   Ele_ESSigRR = 0;
   Ele_SCRawE = 0;
   Ele_SC_ESEnByRawE = 0;
   Ele_HadOverEm = 0;
   Ele_Gen_Pt = 0;
   Ele_Gen_Eta = 0;
   Ele_Gen_Phi = 0;
   Ele_Gen_E = 0;
// Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("iEtaEle1", &iEtaEle1, &b_iEtaEle1);
   fChain->SetBranchAddress("iPhiEle1", &iPhiEle1, &b_iPhiEle1);
   fChain->SetBranchAddress("Hit_ES_Eta_Ele1", &Hit_ES_Eta_Ele1, &b_Hit_ES_Eta_Ele1);
   fChain->SetBranchAddress("Hit_ES_Phi_Ele1", &Hit_ES_Phi_Ele1, &b_Hit_ES_Phi_Ele1);
   fChain->SetBranchAddress("Hit_ES_X_Ele1", &Hit_ES_X_Ele1, &b_Hit_ES_X_Ele1);
   fChain->SetBranchAddress("Hit_ES_Y_Ele1", &Hit_ES_Y_Ele1, &b_Hit_ES_Y_Ele1);
   fChain->SetBranchAddress("Hit_ES_Z_Ele1", &Hit_ES_Z_Ele1, &b_Hit_ES_Z_Ele1);
   fChain->SetBranchAddress("ES_RecHitEnEle1", &ES_RecHitEnEle1, &b_ES_RecHitEnEle1);
   fChain->SetBranchAddress("Hit_Eta_Ele1", &Hit_Eta_Ele1, &b_Hit_Eta_Ele1);
   fChain->SetBranchAddress("Hit_Phi_Ele1", &Hit_Phi_Ele1, &b_Hit_Phi_Ele1);
   fChain->SetBranchAddress("Hit_X_Ele1", &Hit_X_Ele1, &b_Hit_X_Ele1);
   fChain->SetBranchAddress("Hit_Y_Ele1", &Hit_Y_Ele1, &b_Hit_Y_Ele1);
   fChain->SetBranchAddress("Hit_Z_Ele1", &Hit_Z_Ele1, &b_Hit_Z_Ele1);
   fChain->SetBranchAddress("RecHitEnEle1", &RecHitEnEle1, &b_RecHitEnEle1);
   fChain->SetBranchAddress("RecHitFracPho1", &RecHitFracPho1, &b_RecHitFracPho1);
   fChain->SetBranchAddress("RecHitGain1", &RecHitGain1, &b_RecHitGain1);
   fChain->SetBranchAddress("RecHitQuality1", &RecHitQuality1, &b_RecHitQuality1);
   fChain->SetBranchAddress("HitNoiseEle1", &HitNoiseEle1, &b_HitNoiseEle1);
   fChain->SetBranchAddress("iEtaEle2", &iEtaEle2, &b_iEtaEle2);
   fChain->SetBranchAddress("iPhiEle2", &iPhiEle2, &b_iPhiEle2);
   fChain->SetBranchAddress("Hit_ES_Eta_Ele2", &Hit_ES_Eta_Ele2, &b_Hit_ES_Eta_Ele2);
   fChain->SetBranchAddress("Hit_ES_Phi_Ele2", &Hit_ES_Phi_Ele2, &b_Hit_ES_Phi_Ele2);
   fChain->SetBranchAddress("Hit_ES_X_Ele2", &Hit_ES_X_Ele2, &b_Hit_ES_X_Ele2);
   fChain->SetBranchAddress("Hit_ES_Y_Ele2", &Hit_ES_Y_Ele2, &b_Hit_ES_Y_Ele2);
   fChain->SetBranchAddress("Hit_ES_Z_Ele2", &Hit_ES_Z_Ele2, &b_Hit_ES_Z_Ele2);
   fChain->SetBranchAddress("ES_RecHitEnEle2", &ES_RecHitEnEle2, &b_ES_RecHitEnEle2);
   fChain->SetBranchAddress("Hit_Eta_Ele2", &Hit_Eta_Ele2, &b_Hit_Eta_Ele2);
   fChain->SetBranchAddress("Hit_Phi_Ele2", &Hit_Phi_Ele2, &b_Hit_Phi_Ele2);
   fChain->SetBranchAddress("Hit_X_Ele2", &Hit_X_Ele2, &b_Hit_X_Ele2);
   fChain->SetBranchAddress("Hit_Y_Ele2", &Hit_Y_Ele2, &b_Hit_Y_Ele2);
   fChain->SetBranchAddress("Hit_Z_Ele2", &Hit_Z_Ele2, &b_Hit_Z_Ele2);
   fChain->SetBranchAddress("RecHitEnEle2", &RecHitEnEle2, &b_RecHitEnEle2);
   fChain->SetBranchAddress("RecHitFracPho2", &RecHitFracPho2, &b_RecHitFracPho2);
   fChain->SetBranchAddress("RecHitGain2", &RecHitGain2, &b_RecHitGain2);
   fChain->SetBranchAddress("RecHitQuality2", &RecHitQuality2, &b_RecHitQuality2);
   fChain->SetBranchAddress("HitNoiseEle2", &HitNoiseEle2, &b_HitNoiseEle2);
   fChain->SetBranchAddress("iEtaEle3", &iEtaEle3, &b_iEtaEle3);
   fChain->SetBranchAddress("iPhiEle3", &iPhiEle3, &b_iPhiEle3);
   fChain->SetBranchAddress("Hit_ES_Eta_Ele3", &Hit_ES_Eta_Ele3, &b_Hit_ES_Eta_Ele3);
   fChain->SetBranchAddress("Hit_ES_Phi_Ele3", &Hit_ES_Phi_Ele3, &b_Hit_ES_Phi_Ele3);
   fChain->SetBranchAddress("Hit_ES_X_Ele3", &Hit_ES_X_Ele3, &b_Hit_ES_X_Ele3);
   fChain->SetBranchAddress("Hit_ES_Y_Ele3", &Hit_ES_Y_Ele3, &b_Hit_ES_Y_Ele3);
   fChain->SetBranchAddress("Hit_ES_Z_Ele3", &Hit_ES_Z_Ele3, &b_Hit_ES_Z_Ele3);
   fChain->SetBranchAddress("ES_RecHitEnEle3", &ES_RecHitEnEle3, &b_ES_RecHitEnEle3);
   fChain->SetBranchAddress("Hit_Eta_Ele3", &Hit_Eta_Ele3, &b_Hit_Eta_Ele3);
   fChain->SetBranchAddress("Hit_Phi_Ele3", &Hit_Phi_Ele3, &b_Hit_Phi_Ele3);
   fChain->SetBranchAddress("Hit_X_Ele3", &Hit_X_Ele3, &b_Hit_X_Ele3);
   fChain->SetBranchAddress("Hit_Y_Ele3", &Hit_Y_Ele3, &b_Hit_Y_Ele3);
   fChain->SetBranchAddress("Hit_Z_Ele3", &Hit_Z_Ele3, &b_Hit_Z_Ele3);
   fChain->SetBranchAddress("RecHitEnEle3", &RecHitEnEle3, &b_RecHitEnEle3);
   fChain->SetBranchAddress("RecHitFracPho3", &RecHitFracPho3, &b_RecHitFracPho3);
   fChain->SetBranchAddress("RecHitGain3", &RecHitGain3, &b_RecHitGain3);
   fChain->SetBranchAddress("RecHitQuality3", &RecHitQuality3, &b_RecHitQuality3);
   fChain->SetBranchAddress("HitNoiseEle3", &HitNoiseEle3, &b_HitNoiseEle3);
   fChain->SetBranchAddress("iEtaEle4", &iEtaEle4, &b_iEtaEle4);
   fChain->SetBranchAddress("iPhiEle4", &iPhiEle4, &b_iPhiEle4);
   fChain->SetBranchAddress("Hit_ES_Eta_Ele4", &Hit_ES_Eta_Ele4, &b_Hit_ES_Eta_Ele4);
   fChain->SetBranchAddress("Hit_ES_Phi_Ele4", &Hit_ES_Phi_Ele4, &b_Hit_ES_Phi_Ele4);
   fChain->SetBranchAddress("Hit_ES_X_Ele4", &Hit_ES_X_Ele4, &b_Hit_ES_X_Ele4);
   fChain->SetBranchAddress("Hit_ES_Y_Ele4", &Hit_ES_Y_Ele4, &b_Hit_ES_Y_Ele4);
   fChain->SetBranchAddress("Hit_ES_Z_Ele4", &Hit_ES_Z_Ele4, &b_Hit_ES_Z_Ele4);
   fChain->SetBranchAddress("ES_RecHitEnEle4", &ES_RecHitEnEle4, &b_ES_RecHitEnEle4);
   fChain->SetBranchAddress("Hit_Eta_Ele4", &Hit_Eta_Ele4, &b_Hit_Eta_Ele4);
   fChain->SetBranchAddress("Hit_Phi_Ele4", &Hit_Phi_Ele4, &b_Hit_Phi_Ele4);
   fChain->SetBranchAddress("Hit_X_Ele4", &Hit_X_Ele4, &b_Hit_X_Ele4);
   fChain->SetBranchAddress("Hit_Y_Ele4", &Hit_Y_Ele4, &b_Hit_Y_Ele4);
   fChain->SetBranchAddress("Hit_Z_Ele4", &Hit_Z_Ele4, &b_Hit_Z_Ele4);
   fChain->SetBranchAddress("RecHitEnEle4", &RecHitEnEle4, &b_RecHitEnEle4);
   fChain->SetBranchAddress("RecHitFracPho4", &RecHitFracPho4, &b_RecHitFracPho4);
   fChain->SetBranchAddress("RecHitGain4", &RecHitGain4, &b_RecHitGain4);
   fChain->SetBranchAddress("RecHitQuality4", &RecHitQuality4, &b_RecHitQuality4);
   fChain->SetBranchAddress("HitNoiseEle4", &HitNoiseEle4, &b_HitNoiseEle4);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nEle);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("energy_ecal_mustache", &energy_ecal_mustache, &b_energy_ecal_mustache);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("passTightId", &passTightId, &b_passTightId);
   fChain->SetBranchAddress("passMVAMediumId", &passMVAMediumId, &b_passMVAMediumId);
   fChain->SetBranchAddress("Ele_R9", &Ele_R9, &b_Ele_R9);
   fChain->SetBranchAddress("Ele_S4", &Ele_S4, &b_Ele_S4);
   fChain->SetBranchAddress("Ele_SigIEIE", &Ele_SigIEIE, &b_Ele_SigIEIE);
   fChain->SetBranchAddress("Ele_SigIPhiIPhi", &Ele_SigIPhiIPhi, &b_Ele_SigIPhiIPhi);
   fChain->SetBranchAddress("Ele_SCEtaW", &Ele_SCEtaW, &b_Ele_SCEtaW);
   fChain->SetBranchAddress("Ele_SCPhiW", &Ele_SCPhiW, &b_Ele_SCPhiW);
   fChain->SetBranchAddress("Ele_CovIEtaIEta", &Ele_CovIEtaIEta, &b_Ele_CovIEtaIEta);
   fChain->SetBranchAddress("Ele_CovIEtaIPhi", &Ele_CovIEtaIPhi, &b_Ele_CovIEtaIPhi);
   fChain->SetBranchAddress("Ele_ESSigRR", &Ele_ESSigRR, &b_Ele_ESSigRR);
   fChain->SetBranchAddress("Ele_SCRawE", &Ele_SCRawE, &b_Ele_SCRawE);
   fChain->SetBranchAddress("Ele_SC_ESEnByRawE", &Ele_SC_ESEnByRawE, &b_Ele_SC_ESEnByRawE);
   fChain->SetBranchAddress("Ele_HadOverEm", &Ele_HadOverEm, &b_Ele_HadOverEm);
   fChain->SetBranchAddress("Ele_Gen_Pt", &Ele_Gen_Pt, &b_Ele_Gen_Pt);
   fChain->SetBranchAddress("Ele_Gen_Eta", &Ele_Gen_Eta, &b_Ele_Gen_Eta);
   fChain->SetBranchAddress("Ele_Gen_Phi", &Ele_Gen_Phi, &b_Ele_Gen_Phi);
   fChain->SetBranchAddress("Ele_Gen_E", &Ele_Gen_E, &b_Ele_Gen_E);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
  if (!tree)
    return;
  Notify();
  return;
  // End of modification by Somanko
  //  Set branch addresses and branch pointers
  if (!tree2)
    return;
  fChain2 = tree2;
  fCurrent2 = -1;
  fChain2->SetMakeClass(1);

  Notify();
}

Bool_t HGCNtupleVariables::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef HGCNtupleVariables_cxx