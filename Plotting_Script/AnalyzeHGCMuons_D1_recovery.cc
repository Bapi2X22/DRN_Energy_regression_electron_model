#define ANALYZEHGCMuons_D1_cxx

#include "AnalyzeHGCMuons_D1.h"
//#include "TLorentzVector.h"
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;
int ctr = 0;
int main(int argc, char *argv[]) {

  if (argc < 2) {
    cerr << "Please give 3 arguments "
         << "runList "
         << " "
         << "outputFileName"
         << " "
         << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  const char *data = argv[3];
  const char *massP = argv[4];

  cout<<massP<<endl;
  AnalyzeHGCMuons_D1 hgcmuons(inputFileList, outFileName, data,massP);
  // cout << "dataset " << data << " " << endl;

  hgcmuons.EventLoop(data);
  // DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
  // BDT_ratio->Scale(1.0 / BDT_ratio->Integral());

  return 0;
}

double distanceToOne(double ratio) {
    return fabs(ratio - 1.0);
  
}

void AnalyzeHGCMuons_D1::EventLoop(const char *data) {
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
 
  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  int decade = 0;
  int veto_count =0;
  int drnCloserCount = 0;
  int bdtCloserCount = 0;

  for (Long64_t jentry = 0; jentry < nentries; jentry++) {

  
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int(progress);
    if (k > decade)
      // cout << 10 * k << " %" << endl;
      decade = k;

    // ===============read this entry == == == == == == == == == == ==
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

        //Defining gen quantities
    DRN_ratio->SetLineColor(kRed);
    BDT_ratio->SetLineColor(kBlue);
    
    size_t numElectrons = matchedGenEnergy_->size();
      // Flags for barrel and endcap
    bool barrel = true;   // Set to true for barrel condition
    bool endcap = false;  // Set to true for endcap condition
    
    if (numElectrons < 4) continue; // Skip events with fewer than 2 electrons

    // bool skipEvent = false;  // Flag to skip event if the condition is not met

    // if (barrel) {
    //     for (size_t i = 0; i < numElectrons; ++i) {
    //         double etaValue = fabs(eta->at(i));  // Use absolute value of eta
    //         if (etaValue > 1.4) {  // Skip if eta is outside the barrel range
    //             skipEvent = true;
    //             break;
    //         }
    //     }
    // }
    
    // // Apply endcap condition
    // if (endcap) {
    //     for (size_t i = 0; i < numElectrons; ++i) {
    //         double etaValue = fabs(eta->at(i));  // Use absolute value of eta
    //         if (etaValue < 1.4 || etaValue > 2.5) {  // Skip if outside endcap range
    //             skipEvent = true;
    //             break;
    //         }
    //     }
    // }
        
    // // Skip the event if it fails the conditions
    // if (skipEvent) continue;
    
    for (size_t i = 0; i < numElectrons; ++i) {
        double etaValue = fabs(eta->at(i));  // Use absolute value of eta
    
    // Apply barrel or endcap condition to individual electrons
        bool passesEta = false;
    
        if (barrel && etaValue <= 1) {
            passesEta = true;
        }
    
        if (endcap && etaValue >= 1.4 && etaValue <= 2.5) {
            passesEta = true;
        }
    
        // Skip this electron if it doesn't pass the eta condition
        if (!passesEta) continue;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double recoDRN = recoDRNEnergy->at(i);
        double genD = matchedGenEnergy_->at(i);
        double ratioD = recoDRN / genD;
        // cout << "DRN: " << ratioD << endl;
        DRN_ratio->Fill(ratioD);

        double recoBDT = energy_ecal_mustache->at(i);
        double genB = matchedGenEnergy_->at(i);
        double ratioB = recoBDT / genB;
        // cout << "BDT: " << ratioB << endl;
        BDT_ratio->Fill(ratioB);

        // Compare which ratio is closer to 1
        if (distanceToOne(ratioD) < distanceToOne(ratioB)) {
            drnCloserCount++;
        } else {
            bdtCloserCount++;
        }
    }

	}

  // DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
  // BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
  cout << "Total times DRN ratio is closer to 1: " << drnCloserCount << endl;
  cout << "Total times BDT ratio is closer to 1: " << bdtCloserCount << endl;
}
