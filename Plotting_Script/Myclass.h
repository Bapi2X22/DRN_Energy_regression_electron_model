//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 22 10:59:13 2025 by ROOT version 6.32.04
// from TTree T/MyTuple
// found on file: ElectronRecHits_ntuple_3.root
//////////////////////////////////////////////////////////

#ifndef Myclass_h
#define Myclass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class Myclass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxmatchedGenEta = 1;
   static constexpr Int_t kMaxmatchedGenphi = 1;
   static constexpr Int_t kMaxmatchedGenpt = 1;
   static constexpr Int_t kMaxmatchedGenEnergy = 1;

   // Declaration of leaf types
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
   vector<float>   *RecHitFracEle1;
   vector<int>     *RecHitGain1;
   vector<bool>    *RecHitQuality1;
   vector<bool>    *RecHitFlag_kGood_ele1;
   vector<bool>    *RecHitFlag_kPoorReco_ele1;
   vector<bool>    *RecHitFlag_kOutOfTime_ele1;
   vector<bool>    *RecHitFlag_kFaultyHardware_ele1;
   vector<bool>    *RecHitFlag_kNoisy_ele1;
   vector<bool>    *RecHitFlag_kPoorCalib_ele1;
   vector<bool>    *RecHitFlag_kSaturated_ele1;
   vector<bool>    *RecHitFlag_kLeadingEdgeRecovered_ele1;
   vector<bool>    *RecHitFlag_kNeighboursRecovered_ele1;
   vector<bool>    *RecHitFlag_kTowerRecovered_ele1;
   vector<bool>    *RecHitFlag_kDead_ele1;
   vector<bool>    *RecHitFlag_kKilled_ele1;
   vector<bool>    *RecHitFlag_kTPSaturated_ele1;
   vector<bool>    *RecHitFlag_kL1SpikeFlag_ele1;
   vector<bool>    *RecHitFlag_kWeird_ele1;
   vector<bool>    *RecHitFlag_kDiWeird_ele1;
   vector<bool>    *RecHitFlag_kHasSwitchToGain6_ele1;
   vector<bool>    *RecHitFlag_kHasSwitchToGain1_ele1;
   vector<bool>    *RecHitFlag_kESGood_ele1;
   vector<bool>    *RecHitFlag_kESDead_ele1;
   vector<bool>    *RecHitFlag_kESHot_ele1;
   vector<bool>    *RecHitFlag_kESPassBX_ele1;
   vector<bool>    *RecHitFlag_kESTwoGoodRatios_ele1;
   vector<bool>    *RecHitFlag_kESBadRatioFor12_ele1;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Upper_ele1;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Lower_ele1;
   vector<bool>    *RecHitFlag_kESTS1Largest_ele1;
   vector<bool>    *RecHitFlag_kESTS3Largest_ele1;
   vector<bool>    *RecHitFlag_kESTS3Negative_ele1;
   vector<bool>    *RecHitFlag_kESSaturated_ele1;
   vector<bool>    *RecHitFlag_kESTS2Saturated_ele1;
   vector<bool>    *RecHitFlag_kESTS3Saturated_ele1;
   vector<bool>    *RecHitFlag_kESTS13Sigmas_ele1;
   vector<bool>    *RecHitFlag_kESTS15Sigmas_ele1;
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
   vector<float>   *RecHitFracEle2;
   vector<int>     *RecHitGain2;
   vector<bool>    *RecHitQuality2;
   vector<bool>    *RecHitFlag_kGood_ele2;
   vector<bool>    *RecHitFlag_kPoorReco_ele2;
   vector<bool>    *RecHitFlag_kOutOfTime_ele2;
   vector<bool>    *RecHitFlag_kFaultyHardware_ele2;
   vector<bool>    *RecHitFlag_kNoisy_ele2;
   vector<bool>    *RecHitFlag_kPoorCalib_ele2;
   vector<bool>    *RecHitFlag_kSaturated_ele2;
   vector<bool>    *RecHitFlag_kLeadingEdgeRecovered_ele2;
   vector<bool>    *RecHitFlag_kNeighboursRecovered_ele2;
   vector<bool>    *RecHitFlag_kTowerRecovered_ele2;
   vector<bool>    *RecHitFlag_kDead_ele2;
   vector<bool>    *RecHitFlag_kKilled_ele2;
   vector<bool>    *RecHitFlag_kTPSaturated_ele2;
   vector<bool>    *RecHitFlag_kL1SpikeFlag_ele2;
   vector<bool>    *RecHitFlag_kWeird_ele2;
   vector<bool>    *RecHitFlag_kDiWeird_ele2;
   vector<bool>    *RecHitFlag_kHasSwitchToGain6_ele2;
   vector<bool>    *RecHitFlag_kHasSwitchToGain1_ele2;
   vector<bool>    *RecHitFlag_kESGood_ele2;
   vector<bool>    *RecHitFlag_kESDead_ele2;
   vector<bool>    *RecHitFlag_kESHot_ele2;
   vector<bool>    *RecHitFlag_kESPassBX_ele2;
   vector<bool>    *RecHitFlag_kESTwoGoodRatios_ele2;
   vector<bool>    *RecHitFlag_kESBadRatioFor12_ele2;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Upper_ele2;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Lower_ele2;
   vector<bool>    *RecHitFlag_kESTS1Largest_ele2;
   vector<bool>    *RecHitFlag_kESTS3Largest_ele2;
   vector<bool>    *RecHitFlag_kESTS3Negative_ele2;
   vector<bool>    *RecHitFlag_kESSaturated_ele2;
   vector<bool>    *RecHitFlag_kESTS2Saturated_ele2;
   vector<bool>    *RecHitFlag_kESTS3Saturated_ele2;
   vector<bool>    *RecHitFlag_kESTS13Sigmas_ele2;
   vector<bool>    *RecHitFlag_kESTS15Sigmas_ele2;
   Int_t           nElectrons;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *energy;
   vector<float>   *energy_error;
   vector<float>   *energy_ecal_mustache;
   vector<int>     *passLooseId;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
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
   vector<float>   *sumChargedHadronPt;
   vector<float>   *sumChargedParticlePt;
   vector<float>   *sumEcalClusterEt;
   vector<float>   *sumHcalClusterEt;
   vector<float>   *sumNeutralHadronEt;
   vector<float>   *sumPhotonEt;
   vector<float>   *sumPUPt;
   vector<float>   *Ele_EcalPFClusterIso;
   vector<float>   *Ele_HcalPFClusterIso;
   vector<float>   *matchedGenEta_;
   vector<float>   *matchedGenphi_;
   vector<float>   *matchedGenpt_;
   vector<float>   *matchedGenEnergy_;
   vector<float>   *recoDRNEnergy;
   Float_t         rho;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Bool_t          isRefinedSC;

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
   TBranch        *b_RecHitFracEle1;   //!
   TBranch        *b_RecHitGain1;   //!
   TBranch        *b_RecHitQuality1;   //!
   TBranch        *b_RecHitFlag_kGood_ele1;   //!
   TBranch        *b_RecHitFlag_kPoorReco_ele1;   //!
   TBranch        *b_RecHitFlag_kOutOfTime_ele1;   //!
   TBranch        *b_RecHitFlag_kFaultyHardware_ele1;   //!
   TBranch        *b_RecHitFlag_kNoisy_ele1;   //!
   TBranch        *b_RecHitFlag_kPoorCalib_ele1;   //!
   TBranch        *b_RecHitFlag_kSaturated_ele1;   //!
   TBranch        *b_RecHitFlag_kLeadingEdgeRecovered_ele1;   //!
   TBranch        *b_RecHitFlag_kNeighboursRecovered_ele1;   //!
   TBranch        *b_RecHitFlag_kTowerRecovered_ele1;   //!
   TBranch        *b_RecHitFlag_kDead_ele1;   //!
   TBranch        *b_RecHitFlag_kKilled_ele1;   //!
   TBranch        *b_RecHitFlag_kTPSaturated_ele1;   //!
   TBranch        *b_RecHitFlag_kL1SpikeFlag_ele1;   //!
   TBranch        *b_RecHitFlag_kWeird_ele1;   //!
   TBranch        *b_RecHitFlag_kDiWeird_ele1;   //!
   TBranch        *b_RecHitFlag_kHasSwitchToGain6_ele1;   //!
   TBranch        *b_RecHitFlag_kHasSwitchToGain1_ele1;   //!
   TBranch        *b_RecHitFlag_kESGood_ele1;   //!
   TBranch        *b_RecHitFlag_kESDead_ele1;   //!
   TBranch        *b_RecHitFlag_kESHot_ele1;   //!
   TBranch        *b_RecHitFlag_kESPassBX_ele1;   //!
   TBranch        *b_RecHitFlag_kESTwoGoodRatios_ele1;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor12_ele1;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor23Upper_ele1;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor23Lower_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS1Largest_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS3Largest_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS3Negative_ele1;   //!
   TBranch        *b_RecHitFlag_kESSaturated_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS2Saturated_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS3Saturated_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS13Sigmas_ele1;   //!
   TBranch        *b_RecHitFlag_kESTS15Sigmas_ele1;   //!
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
   TBranch        *b_RecHitFracEle2;   //!
   TBranch        *b_RecHitGain2;   //!
   TBranch        *b_RecHitQuality2;   //!
   TBranch        *b_RecHitFlag_kGood_ele2;   //!
   TBranch        *b_RecHitFlag_kPoorReco_ele2;   //!
   TBranch        *b_RecHitFlag_kOutOfTime_ele2;   //!
   TBranch        *b_RecHitFlag_kFaultyHardware_ele2;   //!
   TBranch        *b_RecHitFlag_kNoisy_ele2;   //!
   TBranch        *b_RecHitFlag_kPoorCalib_ele2;   //!
   TBranch        *b_RecHitFlag_kSaturated_ele2;   //!
   TBranch        *b_RecHitFlag_kLeadingEdgeRecovered_ele2;   //!
   TBranch        *b_RecHitFlag_kNeighboursRecovered_ele2;   //!
   TBranch        *b_RecHitFlag_kTowerRecovered_ele2;   //!
   TBranch        *b_RecHitFlag_kDead_ele2;   //!
   TBranch        *b_RecHitFlag_kKilled_ele2;   //!
   TBranch        *b_RecHitFlag_kTPSaturated_ele2;   //!
   TBranch        *b_RecHitFlag_kL1SpikeFlag_ele2;   //!
   TBranch        *b_RecHitFlag_kWeird_ele2;   //!
   TBranch        *b_RecHitFlag_kDiWeird_ele2;   //!
   TBranch        *b_RecHitFlag_kHasSwitchToGain6_ele2;   //!
   TBranch        *b_RecHitFlag_kHasSwitchToGain1_ele2;   //!
   TBranch        *b_RecHitFlag_kESGood_ele2;   //!
   TBranch        *b_RecHitFlag_kESDead_ele2;   //!
   TBranch        *b_RecHitFlag_kESHot_ele2;   //!
   TBranch        *b_RecHitFlag_kESPassBX_ele2;   //!
   TBranch        *b_RecHitFlag_kESTwoGoodRatios_ele2;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor12_ele2;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor23Upper_ele2;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor23Lower_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS1Largest_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS3Largest_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS3Negative_ele2;   //!
   TBranch        *b_RecHitFlag_kESSaturated_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS2Saturated_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS3Saturated_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS13Sigmas_ele2;   //!
   TBranch        *b_RecHitFlag_kESTS15Sigmas_ele2;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_energy_error;   //!
   TBranch        *b_energy_ecal_mustache;   //!
   TBranch        *b_passLooseId;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_passTightId;   //!
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
   TBranch        *b_sumChargedHadronPt;   //!
   TBranch        *b_sumChargedParticlePt;   //!
   TBranch        *b_sumEcalClusterEt;   //!
   TBranch        *b_sumHcalClusterEt;   //!
   TBranch        *b_sumNeutralHadronEt;   //!
   TBranch        *b_sumPhotonEt;   //!
   TBranch        *b_sumPUPt;   //!
   TBranch        *b_Ele_EcalPFClusterIso;   //!
   TBranch        *b_Ele_HcalPFClusterIso;   //!
   TBranch        *b_matchedGenEta_;   //!
   TBranch        *b_matchedGenphi_;   //!
   TBranch        *b_matchedGenpt_;   //!
   TBranch        *b_matchedGenEnergy_;   //!
   TBranch        *b_recoDRNEnergy;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_isRefinedSC;   //!

   Myclass(TTree *tree=0);
   virtual ~Myclass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Myclass_cxx
Myclass::Myclass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ElectronRecHits_ntuple_3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ElectronRecHits_ntuple_3.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ElectronRecHits_ntuple_3.root:/nTuplelize");
      dir->GetObject("T",tree);

   }
   Init(tree);
}

Myclass::~Myclass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Myclass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Myclass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Myclass::Init(TTree *tree)
{
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
   RecHitFracEle1 = 0;
   RecHitGain1 = 0;
   RecHitQuality1 = 0;
   RecHitFlag_kGood_ele1 = 0;
   RecHitFlag_kPoorReco_ele1 = 0;
   RecHitFlag_kOutOfTime_ele1 = 0;
   RecHitFlag_kFaultyHardware_ele1 = 0;
   RecHitFlag_kNoisy_ele1 = 0;
   RecHitFlag_kPoorCalib_ele1 = 0;
   RecHitFlag_kSaturated_ele1 = 0;
   RecHitFlag_kLeadingEdgeRecovered_ele1 = 0;
   RecHitFlag_kNeighboursRecovered_ele1 = 0;
   RecHitFlag_kTowerRecovered_ele1 = 0;
   RecHitFlag_kDead_ele1 = 0;
   RecHitFlag_kKilled_ele1 = 0;
   RecHitFlag_kTPSaturated_ele1 = 0;
   RecHitFlag_kL1SpikeFlag_ele1 = 0;
   RecHitFlag_kWeird_ele1 = 0;
   RecHitFlag_kDiWeird_ele1 = 0;
   RecHitFlag_kHasSwitchToGain6_ele1 = 0;
   RecHitFlag_kHasSwitchToGain1_ele1 = 0;
   RecHitFlag_kESGood_ele1 = 0;
   RecHitFlag_kESDead_ele1 = 0;
   RecHitFlag_kESHot_ele1 = 0;
   RecHitFlag_kESPassBX_ele1 = 0;
   RecHitFlag_kESTwoGoodRatios_ele1 = 0;
   RecHitFlag_kESBadRatioFor12_ele1 = 0;
   RecHitFlag_kESBadRatioFor23Upper_ele1 = 0;
   RecHitFlag_kESBadRatioFor23Lower_ele1 = 0;
   RecHitFlag_kESTS1Largest_ele1 = 0;
   RecHitFlag_kESTS3Largest_ele1 = 0;
   RecHitFlag_kESTS3Negative_ele1 = 0;
   RecHitFlag_kESSaturated_ele1 = 0;
   RecHitFlag_kESTS2Saturated_ele1 = 0;
   RecHitFlag_kESTS3Saturated_ele1 = 0;
   RecHitFlag_kESTS13Sigmas_ele1 = 0;
   RecHitFlag_kESTS15Sigmas_ele1 = 0;
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
   RecHitFracEle2 = 0;
   RecHitGain2 = 0;
   RecHitQuality2 = 0;
   RecHitFlag_kGood_ele2 = 0;
   RecHitFlag_kPoorReco_ele2 = 0;
   RecHitFlag_kOutOfTime_ele2 = 0;
   RecHitFlag_kFaultyHardware_ele2 = 0;
   RecHitFlag_kNoisy_ele2 = 0;
   RecHitFlag_kPoorCalib_ele2 = 0;
   RecHitFlag_kSaturated_ele2 = 0;
   RecHitFlag_kLeadingEdgeRecovered_ele2 = 0;
   RecHitFlag_kNeighboursRecovered_ele2 = 0;
   RecHitFlag_kTowerRecovered_ele2 = 0;
   RecHitFlag_kDead_ele2 = 0;
   RecHitFlag_kKilled_ele2 = 0;
   RecHitFlag_kTPSaturated_ele2 = 0;
   RecHitFlag_kL1SpikeFlag_ele2 = 0;
   RecHitFlag_kWeird_ele2 = 0;
   RecHitFlag_kDiWeird_ele2 = 0;
   RecHitFlag_kHasSwitchToGain6_ele2 = 0;
   RecHitFlag_kHasSwitchToGain1_ele2 = 0;
   RecHitFlag_kESGood_ele2 = 0;
   RecHitFlag_kESDead_ele2 = 0;
   RecHitFlag_kESHot_ele2 = 0;
   RecHitFlag_kESPassBX_ele2 = 0;
   RecHitFlag_kESTwoGoodRatios_ele2 = 0;
   RecHitFlag_kESBadRatioFor12_ele2 = 0;
   RecHitFlag_kESBadRatioFor23Upper_ele2 = 0;
   RecHitFlag_kESBadRatioFor23Lower_ele2 = 0;
   RecHitFlag_kESTS1Largest_ele2 = 0;
   RecHitFlag_kESTS3Largest_ele2 = 0;
   RecHitFlag_kESTS3Negative_ele2 = 0;
   RecHitFlag_kESSaturated_ele2 = 0;
   RecHitFlag_kESTS2Saturated_ele2 = 0;
   RecHitFlag_kESTS3Saturated_ele2 = 0;
   RecHitFlag_kESTS13Sigmas_ele2 = 0;
   RecHitFlag_kESTS15Sigmas_ele2 = 0;
   pt = 0;
   eta = 0;
   phi = 0;
   energy = 0;
   energy_error = 0;
   energy_ecal_mustache = 0;
   passLooseId = 0;
   passMediumId = 0;
   passTightId = 0;
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
   sumChargedHadronPt = 0;
   sumChargedParticlePt = 0;
   sumEcalClusterEt = 0;
   sumHcalClusterEt = 0;
   sumNeutralHadronEt = 0;
   sumPhotonEt = 0;
   sumPUPt = 0;
   Ele_EcalPFClusterIso = 0;
   Ele_HcalPFClusterIso = 0;
   matchedGenEta_ = 0;
   matchedGenphi_ = 0;
   matchedGenpt_ = 0;
   matchedGenEnergy_ = 0;
   recoDRNEnergy = 0;
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
   fChain->SetBranchAddress("RecHitFracEle1", &RecHitFracEle1, &b_RecHitFracEle1);
   fChain->SetBranchAddress("RecHitGain1", &RecHitGain1, &b_RecHitGain1);
   fChain->SetBranchAddress("RecHitQuality1", &RecHitQuality1, &b_RecHitQuality1);
   fChain->SetBranchAddress("RecHitFlag_kGood_ele1", &RecHitFlag_kGood_ele1, &b_RecHitFlag_kGood_ele1);
   fChain->SetBranchAddress("RecHitFlag_kPoorReco_ele1", &RecHitFlag_kPoorReco_ele1, &b_RecHitFlag_kPoorReco_ele1);
   fChain->SetBranchAddress("RecHitFlag_kOutOfTime_ele1", &RecHitFlag_kOutOfTime_ele1, &b_RecHitFlag_kOutOfTime_ele1);
   fChain->SetBranchAddress("RecHitFlag_kFaultyHardware_ele1", &RecHitFlag_kFaultyHardware_ele1, &b_RecHitFlag_kFaultyHardware_ele1);
   fChain->SetBranchAddress("RecHitFlag_kNoisy_ele1", &RecHitFlag_kNoisy_ele1, &b_RecHitFlag_kNoisy_ele1);
   fChain->SetBranchAddress("RecHitFlag_kPoorCalib_ele1", &RecHitFlag_kPoorCalib_ele1, &b_RecHitFlag_kPoorCalib_ele1);
   fChain->SetBranchAddress("RecHitFlag_kSaturated_ele1", &RecHitFlag_kSaturated_ele1, &b_RecHitFlag_kSaturated_ele1);
   fChain->SetBranchAddress("RecHitFlag_kLeadingEdgeRecovered_ele1", &RecHitFlag_kLeadingEdgeRecovered_ele1, &b_RecHitFlag_kLeadingEdgeRecovered_ele1);
   fChain->SetBranchAddress("RecHitFlag_kNeighboursRecovered_ele1", &RecHitFlag_kNeighboursRecovered_ele1, &b_RecHitFlag_kNeighboursRecovered_ele1);
   fChain->SetBranchAddress("RecHitFlag_kTowerRecovered_ele1", &RecHitFlag_kTowerRecovered_ele1, &b_RecHitFlag_kTowerRecovered_ele1);
   fChain->SetBranchAddress("RecHitFlag_kDead_ele1", &RecHitFlag_kDead_ele1, &b_RecHitFlag_kDead_ele1);
   fChain->SetBranchAddress("RecHitFlag_kKilled_ele1", &RecHitFlag_kKilled_ele1, &b_RecHitFlag_kKilled_ele1);
   fChain->SetBranchAddress("RecHitFlag_kTPSaturated_ele1", &RecHitFlag_kTPSaturated_ele1, &b_RecHitFlag_kTPSaturated_ele1);
   fChain->SetBranchAddress("RecHitFlag_kL1SpikeFlag_ele1", &RecHitFlag_kL1SpikeFlag_ele1, &b_RecHitFlag_kL1SpikeFlag_ele1);
   fChain->SetBranchAddress("RecHitFlag_kWeird_ele1", &RecHitFlag_kWeird_ele1, &b_RecHitFlag_kWeird_ele1);
   fChain->SetBranchAddress("RecHitFlag_kDiWeird_ele1", &RecHitFlag_kDiWeird_ele1, &b_RecHitFlag_kDiWeird_ele1);
   fChain->SetBranchAddress("RecHitFlag_kHasSwitchToGain6_ele1", &RecHitFlag_kHasSwitchToGain6_ele1, &b_RecHitFlag_kHasSwitchToGain6_ele1);
   fChain->SetBranchAddress("RecHitFlag_kHasSwitchToGain1_ele1", &RecHitFlag_kHasSwitchToGain1_ele1, &b_RecHitFlag_kHasSwitchToGain1_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESGood_ele1", &RecHitFlag_kESGood_ele1, &b_RecHitFlag_kESGood_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESDead_ele1", &RecHitFlag_kESDead_ele1, &b_RecHitFlag_kESDead_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESHot_ele1", &RecHitFlag_kESHot_ele1, &b_RecHitFlag_kESHot_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESPassBX_ele1", &RecHitFlag_kESPassBX_ele1, &b_RecHitFlag_kESPassBX_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTwoGoodRatios_ele1", &RecHitFlag_kESTwoGoodRatios_ele1, &b_RecHitFlag_kESTwoGoodRatios_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor12_ele1", &RecHitFlag_kESBadRatioFor12_ele1, &b_RecHitFlag_kESBadRatioFor12_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor23Upper_ele1", &RecHitFlag_kESBadRatioFor23Upper_ele1, &b_RecHitFlag_kESBadRatioFor23Upper_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor23Lower_ele1", &RecHitFlag_kESBadRatioFor23Lower_ele1, &b_RecHitFlag_kESBadRatioFor23Lower_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS1Largest_ele1", &RecHitFlag_kESTS1Largest_ele1, &b_RecHitFlag_kESTS1Largest_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Largest_ele1", &RecHitFlag_kESTS3Largest_ele1, &b_RecHitFlag_kESTS3Largest_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Negative_ele1", &RecHitFlag_kESTS3Negative_ele1, &b_RecHitFlag_kESTS3Negative_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESSaturated_ele1", &RecHitFlag_kESSaturated_ele1, &b_RecHitFlag_kESSaturated_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS2Saturated_ele1", &RecHitFlag_kESTS2Saturated_ele1, &b_RecHitFlag_kESTS2Saturated_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Saturated_ele1", &RecHitFlag_kESTS3Saturated_ele1, &b_RecHitFlag_kESTS3Saturated_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS13Sigmas_ele1", &RecHitFlag_kESTS13Sigmas_ele1, &b_RecHitFlag_kESTS13Sigmas_ele1);
   fChain->SetBranchAddress("RecHitFlag_kESTS15Sigmas_ele1", &RecHitFlag_kESTS15Sigmas_ele1, &b_RecHitFlag_kESTS15Sigmas_ele1);
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
   fChain->SetBranchAddress("RecHitFracEle2", &RecHitFracEle2, &b_RecHitFracEle2);
   fChain->SetBranchAddress("RecHitGain2", &RecHitGain2, &b_RecHitGain2);
   fChain->SetBranchAddress("RecHitQuality2", &RecHitQuality2, &b_RecHitQuality2);
   fChain->SetBranchAddress("RecHitFlag_kGood_ele2", &RecHitFlag_kGood_ele2, &b_RecHitFlag_kGood_ele2);
   fChain->SetBranchAddress("RecHitFlag_kPoorReco_ele2", &RecHitFlag_kPoorReco_ele2, &b_RecHitFlag_kPoorReco_ele2);
   fChain->SetBranchAddress("RecHitFlag_kOutOfTime_ele2", &RecHitFlag_kOutOfTime_ele2, &b_RecHitFlag_kOutOfTime_ele2);
   fChain->SetBranchAddress("RecHitFlag_kFaultyHardware_ele2", &RecHitFlag_kFaultyHardware_ele2, &b_RecHitFlag_kFaultyHardware_ele2);
   fChain->SetBranchAddress("RecHitFlag_kNoisy_ele2", &RecHitFlag_kNoisy_ele2, &b_RecHitFlag_kNoisy_ele2);
   fChain->SetBranchAddress("RecHitFlag_kPoorCalib_ele2", &RecHitFlag_kPoorCalib_ele2, &b_RecHitFlag_kPoorCalib_ele2);
   fChain->SetBranchAddress("RecHitFlag_kSaturated_ele2", &RecHitFlag_kSaturated_ele2, &b_RecHitFlag_kSaturated_ele2);
   fChain->SetBranchAddress("RecHitFlag_kLeadingEdgeRecovered_ele2", &RecHitFlag_kLeadingEdgeRecovered_ele2, &b_RecHitFlag_kLeadingEdgeRecovered_ele2);
   fChain->SetBranchAddress("RecHitFlag_kNeighboursRecovered_ele2", &RecHitFlag_kNeighboursRecovered_ele2, &b_RecHitFlag_kNeighboursRecovered_ele2);
   fChain->SetBranchAddress("RecHitFlag_kTowerRecovered_ele2", &RecHitFlag_kTowerRecovered_ele2, &b_RecHitFlag_kTowerRecovered_ele2);
   fChain->SetBranchAddress("RecHitFlag_kDead_ele2", &RecHitFlag_kDead_ele2, &b_RecHitFlag_kDead_ele2);
   fChain->SetBranchAddress("RecHitFlag_kKilled_ele2", &RecHitFlag_kKilled_ele2, &b_RecHitFlag_kKilled_ele2);
   fChain->SetBranchAddress("RecHitFlag_kTPSaturated_ele2", &RecHitFlag_kTPSaturated_ele2, &b_RecHitFlag_kTPSaturated_ele2);
   fChain->SetBranchAddress("RecHitFlag_kL1SpikeFlag_ele2", &RecHitFlag_kL1SpikeFlag_ele2, &b_RecHitFlag_kL1SpikeFlag_ele2);
   fChain->SetBranchAddress("RecHitFlag_kWeird_ele2", &RecHitFlag_kWeird_ele2, &b_RecHitFlag_kWeird_ele2);
   fChain->SetBranchAddress("RecHitFlag_kDiWeird_ele2", &RecHitFlag_kDiWeird_ele2, &b_RecHitFlag_kDiWeird_ele2);
   fChain->SetBranchAddress("RecHitFlag_kHasSwitchToGain6_ele2", &RecHitFlag_kHasSwitchToGain6_ele2, &b_RecHitFlag_kHasSwitchToGain6_ele2);
   fChain->SetBranchAddress("RecHitFlag_kHasSwitchToGain1_ele2", &RecHitFlag_kHasSwitchToGain1_ele2, &b_RecHitFlag_kHasSwitchToGain1_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESGood_ele2", &RecHitFlag_kESGood_ele2, &b_RecHitFlag_kESGood_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESDead_ele2", &RecHitFlag_kESDead_ele2, &b_RecHitFlag_kESDead_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESHot_ele2", &RecHitFlag_kESHot_ele2, &b_RecHitFlag_kESHot_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESPassBX_ele2", &RecHitFlag_kESPassBX_ele2, &b_RecHitFlag_kESPassBX_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTwoGoodRatios_ele2", &RecHitFlag_kESTwoGoodRatios_ele2, &b_RecHitFlag_kESTwoGoodRatios_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor12_ele2", &RecHitFlag_kESBadRatioFor12_ele2, &b_RecHitFlag_kESBadRatioFor12_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor23Upper_ele2", &RecHitFlag_kESBadRatioFor23Upper_ele2, &b_RecHitFlag_kESBadRatioFor23Upper_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor23Lower_ele2", &RecHitFlag_kESBadRatioFor23Lower_ele2, &b_RecHitFlag_kESBadRatioFor23Lower_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS1Largest_ele2", &RecHitFlag_kESTS1Largest_ele2, &b_RecHitFlag_kESTS1Largest_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Largest_ele2", &RecHitFlag_kESTS3Largest_ele2, &b_RecHitFlag_kESTS3Largest_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Negative_ele2", &RecHitFlag_kESTS3Negative_ele2, &b_RecHitFlag_kESTS3Negative_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESSaturated_ele2", &RecHitFlag_kESSaturated_ele2, &b_RecHitFlag_kESSaturated_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS2Saturated_ele2", &RecHitFlag_kESTS2Saturated_ele2, &b_RecHitFlag_kESTS2Saturated_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Saturated_ele2", &RecHitFlag_kESTS3Saturated_ele2, &b_RecHitFlag_kESTS3Saturated_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS13Sigmas_ele2", &RecHitFlag_kESTS13Sigmas_ele2, &b_RecHitFlag_kESTS13Sigmas_ele2);
   fChain->SetBranchAddress("RecHitFlag_kESTS15Sigmas_ele2", &RecHitFlag_kESTS15Sigmas_ele2, &b_RecHitFlag_kESTS15Sigmas_ele2);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nEle);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("energy_error", &energy_error, &b_energy_error);
   fChain->SetBranchAddress("energy_ecal_mustache", &energy_ecal_mustache, &b_energy_ecal_mustache);
   fChain->SetBranchAddress("passLooseId", &passLooseId, &b_passLooseId);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("passTightId", &passTightId, &b_passTightId);
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
   fChain->SetBranchAddress("sumChargedHadronPt", &sumChargedHadronPt, &b_sumChargedHadronPt);
   fChain->SetBranchAddress("sumChargedParticlePt", &sumChargedParticlePt, &b_sumChargedParticlePt);
   fChain->SetBranchAddress("sumEcalClusterEt", &sumEcalClusterEt, &b_sumEcalClusterEt);
   fChain->SetBranchAddress("sumHcalClusterEt", &sumHcalClusterEt, &b_sumHcalClusterEt);
   fChain->SetBranchAddress("sumNeutralHadronEt", &sumNeutralHadronEt, &b_sumNeutralHadronEt);
   fChain->SetBranchAddress("sumPhotonEt", &sumPhotonEt, &b_sumPhotonEt);
   fChain->SetBranchAddress("sumPUPt", &sumPUPt, &b_sumPUPt);
   fChain->SetBranchAddress("Ele_EcalPFClusterIso", &Ele_EcalPFClusterIso, &b_Ele_EcalPFClusterIso);
   fChain->SetBranchAddress("Ele_HcalPFClusterIso", &Ele_HcalPFClusterIso, &b_Ele_HcalPFClusterIso);
   fChain->SetBranchAddress("matchedGenEta_", &matchedGenEta_, &b_matchedGenEta_);
   fChain->SetBranchAddress("matchedGenphi_", &matchedGenphi_, &b_matchedGenphi_);
   fChain->SetBranchAddress("matchedGenpt_", &matchedGenpt_, &b_matchedGenpt_);
   fChain->SetBranchAddress("matchedGenEnergy_", &matchedGenEnergy_, &b_matchedGenEnergy_);
   fChain->SetBranchAddress("recoDRNEnergy", &recoDRNEnergy, &b_recoDRNEnergy);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("isRefinedSC", &isRefinedSC, &b_isRefinedSC);
   Notify();
}

bool Myclass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void Myclass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Myclass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Myclass_cxx
