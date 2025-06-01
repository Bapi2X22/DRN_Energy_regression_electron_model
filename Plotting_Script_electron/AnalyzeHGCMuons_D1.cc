#define ANALYZEHGCMuons_D1_cxx

#include "AnalyzeHGCMuons_D1.h"
//#include "TLorentzVector.h"
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream> // Required for file operations

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

  std::ofstream outFile("lorentz_vectors.txt", std::ios::app); // appends to the file
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
    
   // size_t numElectrons = matchedGenEnergy_->size();
    size_t numElectrons = pt->size();
    if (numElectrons!=2) continue;

    bool pT = true;
    bool BoolR9 = true;
    bool pTcut =true;
      // Flags for barrel and endcap
    bool barrel = false;   // Set to true for barrel condition
    bool endcap = true;  // Set to true for endcap condition

    // bool eta1 = false
    // bool eta2 = false
    // bool eta3 = false
    // bool eta4 = false
    // bool eta5 = false

      
    // if (numElectrons < 4) continue; // Skip events with fewer than 2 electrons

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
    double pT_val1;
    double eta_val1;
    double phi_val1;
    double pT_val2;
    double eta_val2;
    double phi_val2;
    double Energy1_DRN;
    double Energy1_BDT;
    double Energy2_DRN;
    double Energy2_BDT;
    
    int iter = 0; 
    for (size_t i = 0; i < numElectrons; ++i) {
        if (iter == 2){
            continue;
        };
        double etaValue = fabs(eta->at(i));  // Use absolute value of eta
        double pTValue = pt->at(i);
    
    // Apply barrel or endcap condition to individual electrons
        bool passesEta = false;
        bool passpT = false;
    
        if (barrel && etaValue <= 1.442) {
            passesEta = true;
        }
    
        if (endcap && etaValue >= 1.566 && etaValue <= 2.5) {
            passesEta = true;
        }

        if (pTcut && pTValue >= 20.0 && pTValue <= 300.0) {
            passpT = true;
        }
    
        // Skip this electron if it doesn't pass the eta condition
        if (!passesEta) continue; // [Uncomment this if you need to apply eta condition]
        if (!passpT) continue;


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


      //  double pT_val1;
      //  double eta_val1;
      //  double phi_val1;
      //  double pT_val2;
      //  double eta_val2;
      //  double phi_val2;
      //  double Energy1_DRN;
      //  double Energy1_BDT;
      //  double Energy2_DRN;
      //  double Energy2_BDT;


///////////////////////////////////getting the number of Rechits/////////////////////////////////////////////
        // size_t numRec1 = RecHitEnEle1->size();
        // size_t numRec2 = RecHitEnEle2->size();
        // size_t numRec3 = RecHitEnEle3->size();
        // size_t numRec4 = RecHitEnEle4->size();

///////////////////////////////////finding the invarient mass ditribution////////////////////////////////////
        if (numElectrons>=2) {
            if (i==0){
                pT_val1 = pt->at(i);
                eta_val1 = eta->at(i);
                phi_val1 = phi->at(i);
                Energy1_DRN = recoDRNEnergy->at(i);
                Energy1_BDT = energy_ecal_mustache->at(i);
            }
            if (i==1){
                pT_val2 = pt->at(i);
                eta_val2 = eta->at(i);
                phi_val2 = phi->at(i);
                Energy2_DRN = recoDRNEnergy->at(i);
                Energy2_BDT = energy_ecal_mustache->at(i);
            }
        }
        

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
            size_t numRec1 = RecHitEnEle1->size();
            Num_rechits->Fill(numRec1);
            double pT = pt->at(i);
            if (pT){
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }
            double R9_values = Ele_R9->at(i);
            if (BoolR9){
            for (size_t k = 0; k < R9_bins.size(); ++k) {
                if (R9_values >= R9_bins[k].first && R9_values < R9_bins[k].second) {
                    R9_histos_BDT[k]->Fill(ratioB);
                    R9_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }
            double etaValues = fabs(eta->at(i));
            for (size_t k = 0; k < EE_eta_bins.size(); ++k) {
                if (etaValues >= EE_eta_bins[k].first && etaValues < EE_eta_bins[k].second) {
                    eta_histos_BDT_EE[k]->Fill(ratioB);
                    eta_histos_DRN_EE[k]->Fill(ratioD);
                    break;
        }
    }
            for (size_t k = 0; k < EB_eta_bins.size(); ++k) {
                if (etaValues >= EB_eta_bins[k].first && etaValues < EB_eta_bins[k].second) {
                    eta_histos_BDT_EB[k]->Fill(ratioB);
                    eta_histos_DRN_EB[k]->Fill(ratioD);
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
            size_t numRec2 = RecHitEnEle2->size();
            Num_rechits->Fill(numRec2);
            double pT = pt->at(i);
            if (pT){
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }
            double R9_values = Ele_R9->at(i);
            if (BoolR9){
            for (size_t k = 0; k < R9_bins.size(); ++k) {
                if (R9_values >= R9_bins[k].first && R9_values < R9_bins[k].second) {
                    R9_histos_BDT[k]->Fill(ratioB);
                    R9_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }

            double etaValues = fabs(eta->at(i));
            for (size_t k = 0; k < EE_eta_bins.size(); ++k) {
                if (etaValues >= EE_eta_bins[k].first && etaValues < EE_eta_bins[k].second) {
                    eta_histos_BDT_EE[k]->Fill(ratioB);
                    eta_histos_DRN_EE[k]->Fill(ratioD);
                    break;
        }
    }
            for (size_t k = 0; k < EB_eta_bins.size(); ++k) {
                if (etaValues >= EB_eta_bins[k].first && etaValues < EB_eta_bins[k].second) {
                    eta_histos_BDT_EB[k]->Fill(ratioB);
                    eta_histos_DRN_EB[k]->Fill(ratioD);
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
            size_t numRec3 = RecHitEnEle3->size();
            Num_rechits->Fill(numRec3);
            double pT = pt->at(i);
            if (pT){
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }
            double R9_values = Ele_R9->at(i);
            if (BoolR9){
            for (size_t k = 0; k < R9_bins.size(); ++k) {
                if (R9_values >= R9_bins[k].first && R9_values < R9_bins[k].second) {
                    R9_histos_BDT[k]->Fill(ratioB);
                    R9_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }

            double etaValues = fabs(eta->at(i));
            for (size_t k = 0; k < EE_eta_bins.size(); ++k) {
                if (etaValues >= EE_eta_bins[k].first && etaValues < EE_eta_bins[k].second) {
                    eta_histos_BDT_EE[k]->Fill(ratioB);
                    eta_histos_DRN_EE[k]->Fill(ratioD);
                    break;
        }
    }
            for (size_t k = 0; k < EB_eta_bins.size(); ++k) {
                if (etaValues >= EB_eta_bins[k].first && etaValues < EB_eta_bins[k].second) {
                    eta_histos_BDT_EB[k]->Fill(ratioB);
                    eta_histos_DRN_EB[k]->Fill(ratioD);
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
            size_t numRec4 = RecHitEnEle4->size();
            Num_rechits->Fill(numRec4);
            double pT = pt->at(i);
            if (pT){
            for (size_t k = 0; k < pT_bins.size(); ++k) {
                if (pT >= pT_bins[k].first && pT < pT_bins[k].second) {
                    pT_histos_BDT[k]->Fill(ratioB);
                    pT_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }
            double R9_values = Ele_R9->at(i);
            if (BoolR9){
            for (size_t k = 0; k < R9_bins.size(); ++k) {
                if (R9_values >= R9_bins[k].first && R9_values < R9_bins[k].second) {
                    R9_histos_BDT[k]->Fill(ratioB);
                    R9_histos_DRN[k]->Fill(ratioD);
                    break;
        }
    }
            }

            double etaValues = fabs(eta->at(i));
            for (size_t k = 0; k < EE_eta_bins.size(); ++k) {
                if (etaValues >= EE_eta_bins[k].first && etaValues < EE_eta_bins[k].second) {
                    eta_histos_BDT_EE[k]->Fill(ratioB);
                    eta_histos_DRN_EE[k]->Fill(ratioD);
                    break;
        }
    }
            for (size_t k = 0; k < EB_eta_bins.size(); ++k) {
                if (etaValues >= EB_eta_bins[k].first && etaValues < EB_eta_bins[k].second) {
                    eta_histos_BDT_EB[k]->Fill(ratioB);
                    eta_histos_DRN_EB[k]->Fill(ratioD);
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

        // Num_rechits->Fill(numRec1);
        // Num_rechits->Fill(numRec2);
        // Num_rechits->Fill(numRec3);
        // Num_rechits->Fill(numRec4);

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

        iter++;
    } 
      if (numElectrons==2){
        TLorentzVector eplus_BDT, eminus_BDT, eplus_DRN, eminus_DRN;
        eplus_BDT.SetPtEtaPhiE(pT_val1, eta_val1, phi_val1, Energy1_BDT);
        eminus_BDT.SetPtEtaPhiE(pT_val2, eta_val2, phi_val2, Energy2_BDT);

        eplus_DRN.SetPtEtaPhiE(pT_val1, eta_val1, phi_val1, Energy1_DRN);
        eminus_DRN.SetPtEtaPhiE(pT_val2, eta_val2, phi_val2, Energy2_DRN);

        // Sum the two
        TLorentzVector total_BDT = eplus_BDT + eminus_BDT;
        double mass_BDT = total_BDT.M();
        mee_BDT->Fill(mass_BDT);
        TLorentzVector total_DRN = eplus_DRN + eminus_DRN;
        double mass_DRN = total_DRN.M();
        mee_DRN->Fill(mass_DRN);
      
    if (outFile.is_open()) {
        outFile << "BDT: eplus (" 
                << eplus_BDT.Pt() << ", " << eplus_BDT.Eta() << ", " << eplus_BDT.Phi() << ", " << eplus_BDT.M() << "), "
                << "eminus (" 
                << eminus_BDT.Pt() << ", " << eminus_BDT.Eta() << ", " << eminus_BDT.Phi() << ", " << eminus_BDT.M() << "), "
                << "mass = " << mass_BDT << "\n";

        outFile << "DRN: eplus (" 
                << eplus_DRN.Pt() << ", " << eplus_DRN.Eta() << ", " << eplus_DRN.Phi() << ", " << eplus_DRN.M() << "), "
                << "eminus (" 
                << eminus_DRN.Pt() << ", " << eminus_DRN.Eta() << ", " << eminus_DRN.Phi() << ", " << eminus_DRN.M() << "), "
                << "mass = " << mass_DRN << "\n\n";
    }
	}

//      outFile.close();
}
  // DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
  // BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
  cout << "Total times DRN ratio is closer to 1: " << drnCloserCount << endl;
  cout << "Total times BDT ratio is closer to 1: " << bdtCloserCount << endl;
  outFile.close();
}
