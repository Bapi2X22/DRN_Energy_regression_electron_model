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
    bool barrel = false;   // Set to true for barrel condition
    bool endcap = true;  // Set to true for endcap condition
    
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
    
        if (barrel && etaValue <= 1.4) {
            passesEta = true;
        }
    
        if (endcap && etaValue >= 1.5 && etaValue <= 2.4) {
            passesEta = true;
        }
    
        // Skip this electron if it doesn't pass the eta condition
        if (!passesEta) continue; // [Uncomment this if you need to apply eta condition]


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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double R9_values = Ele_R9->at(i);
        double SigIEIE_values = Ele_SigIEIE->at(i);
///////////////////////////////////getting the number of Rechits/////////////////////////////////////////////
        size_t numRec1 = RecHitEnEle1->size();
        size_t numRec2 = RecHitEnEle2->size();
        size_t numRec3 = RecHitEnEle3->size();
        size_t numRec4 = RecHitEnEle4->size();

/////////////////////////////////////Accessing the energy of the rechits//////////////////////////////////////

        // Loop over all rechits in each electron
        if (i==0){
        double totalEnergyEle1 = 0.0;
        for (size_t j = 0; j < RecHitEnEle1->size(); ++j) {
            double energyEle1 = RecHitEnEle1->at(j);
            double HitNoise1 = HitNoiseEle1->at(j);
            totalEnergyEle1 += energyEle1;
            Energy_rechits->Fill(energyEle1);
            HitNoise->Fill(HitNoise1);
        }
            Sum_En_rechits->Fill(totalEnergyEle1);
            double pT = pt->at(i);
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
        }

        if (i==1){
        double totalEnergyEle2 = 0.0;
        for (size_t j = 0; j < RecHitEnEle2->size(); ++j) {
            double energyEle2 = RecHitEnEle2->at(j);
            double HitNoise2 = HitNoiseEle2->at(j);
            totalEnergyEle2 += energyEle2;
            Energy_rechits->Fill(energyEle2);
            HitNoise->Fill(HitNoise2);
        }
            Sum_En_rechits->Fill(totalEnergyEle2);
            double pT = pt->at(i);
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
        }

        if (i==2){
        double totalEnergyEle3 = 0.0;
        for (size_t j = 0; j < RecHitEnEle3->size(); ++j) {
            double energyEle3 = RecHitEnEle3->at(j);
            double HitNoise3 = HitNoiseEle3->at(j);
            totalEnergyEle3 += energyEle3;
            Energy_rechits->Fill(energyEle3);
            HitNoise->Fill(HitNoise3);
        }
            Sum_En_rechits->Fill(totalEnergyEle3);
            double pT = pt->at(i);
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
        }

        if (i==3){
        double totalEnergyEle4 = 0.0;
        for (size_t j = 0; j < RecHitEnEle4->size(); ++j) {
            double energyEle4 = RecHitEnEle4->at(j);
            double HitNoise4 = HitNoiseEle4->at(j);
            totalEnergyEle4 += energyEle4;
            Energy_rechits->Fill(energyEle4);
            HitNoise->Fill(HitNoise4);
        }
            Sum_En_rechits->Fill(totalEnergyEle4);
            double pT = pt->at(i);
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
        }
        

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // double ratioD = recoDRN / genD;
        // // cout << "DRN: " << ratioD << endl;
        // DRN_ratio->Fill(ratioD);

        // double recoBDT = energy_ecal_mustache->at(i);
        // double genB = matchedGenEnergy_->at(i);
        // double ratioB = recoBDT / genB;
        // // cout << "BDT: " << ratioB << endl;
        // BDT_ratio->Fill(ratioB);

        Num_rechits->Fill(numRec1);
        Num_rechits->Fill(numRec2);
        Num_rechits->Fill(numRec3);
        Num_rechits->Fill(numRec4);

        // Sum_En_rechits->Fill(totalEnergyEle1);
        // Sum_En_rechits->Fill(totalEnergyEle2);
        // Sum_En_rechits->Fill(totalEnergyEle3);
        // Sum_En_rechits->Fill(totalEnergyEle4);

        R9->Fill(R9_values);
        SigIEIE->Fill(SigIEIE_values);
    

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
