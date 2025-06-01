#define ANALYZEHGCMuons_cxx

#include "AnalyzeHGCMuons.h"
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
  AnalyzeHGCMuons hgcmuons(inputFileList, outFileName, data,massP);
  // cout << "dataset " << data << " " << endl;

  hgcmuons.EventLoop(data);
  // DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
  // BDT_ratio->Scale(1.0 / BDT_ratio->Integral());

  return 0;
}
void AnalyzeHGCMuons::EventLoop(const char *data) {
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
 
  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  int decade = 0;
   int veto_count =0;

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

//     if(H_Gen_mass->size()==0){continue;}
//    if(A_lead_Pho_Gen_Pt->size() <2){
// 	cout<<"A1 eta"<<"\t"<<A_lead_Gen_eta->at(0)<<endl;
// 	cout<<"A1 eta"<<"\t"<<A_lead_Gen_pt->at(0)<<endl;
// 	veto_count = veto_count+1;
// 		}
//     if(A_lead_Pho_Gen_Pt->size() <2 || A_sublead_Pho_Gen_Pt->size()<2){
// 	continue;} 
//     fillhist = 0;
//         pHgen.SetPtEtaPhiM(H_Gen_pt->at(0), H_Gen_eta->at(0), H_Gen_phi->at(0), H_Gen_mass->at(0));

// 	if(A_lead_Pho_Gen_Pt->size() <2)
// 	  {
// 		continue;
// 	  }
//if (jentry==10){break;}
//cout<< "Entry no." << "\t"<<jentry<< "\t" <<"No. of Gen photons"<<"\t"<< A_lead_Pho_Gen_Pt->size()+A_sublead_Pho_Gen_Pt->size()<<endl;
        //Defining gen quantities
        DRN_ratio->SetLineColor(kRed);
        BDT_ratio->SetLineColor(kBlue);

        // float GenD1 = matchedGenEnergy_->at(0);
    if(recoDRNEnergy->size()==0){continue;}
    // if(matchedGenEnergy_->size()==0){continue;}
    // if(energy_ecal_mustache->size()==0){continue;}
    if(recoDRNEnergy->size()==1){continue;}
    //   float Reco_DRN1 = recoDRNEnergy->at(0);
    //   float GenD1 = matchedGenEnergy_->at(0);
    //   float ratioD1 = Reco_DRN1/GenD1;
    //   // cout << Reco_DRN1 <<endl;
    //   cout<< "DRN: "<<ratioD1<<endl;
    //   DRN_ratio->Fill(ratioD1);
    // }
    
      if(recoDRNEnergy->size()==2){
        float Reco_DRN2 = recoDRNEnergy->at(0);
        float GenD2 = matchedGenEnergy_->at(0);
        float ratioD2 = Reco_DRN2/GenD2;
        cout<< "DRN: "<<ratioD2<<endl;
        DRN_ratio->Fill(ratioD2);
        float Reco_DRN3 = recoDRNEnergy->at(1);
        float GenD3 = matchedGenEnergy_->at(1);
        float ratioD3 = Reco_DRN3/GenD3;
        cout<< "DRN: "<<ratioD3<<endl;
        DRN_ratio->Fill(ratioD3);
        // DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
      }

      if(energy_ecal_mustache->size()==0){continue;}
      // if(matchedGenEnergy_->size()==0){continue;}
      // if(energy_ecal_mustache->size()==0){continue;}
      if(energy_ecal_mustache->size()==1){continue;}
      //   float Reco_BDT1 = energy_ecal_mustache->at(0);
      //   float GenB1 = matchedGenEnergy_->at(0);
      //   float ratioB1 = Reco_BDT1/GenB1;
      //   cout<< "BDT: "<<ratioB1<<endl;
      //   // cout << Reco_DRN1 <<endl;
      //   BDT_ratio->Fill(ratioB1);
      // }
      
        if(energy_ecal_mustache->size()==2){
          float Reco_BDT2 = energy_ecal_mustache ->at(0);
          float GenB2 = matchedGenEnergy_->at(0);
          float ratioB2 = Reco_BDT2/GenB2;
          cout<< "BDT: "<<ratioB2<<endl;
          BDT_ratio->Fill(ratioB2);
          float Reco_BDT3 = energy_ecal_mustache->at(1);
          float GenB3 = matchedGenEnergy_->at(1);
          float ratioB3 = Reco_BDT3/GenB3;
          cout<< "BDT: "<<ratioB3<<endl;
          BDT_ratio->Fill(ratioB3);
          // BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
  
        }


// TCanvas *c1 = new TCanvas("c1", "Reco/Gen Ratios", 800, 600);
// DRN_ratio->Draw("HIST");     // Draw DRN histogram
// BDT_ratio->Draw("HIST SAME"); // Overlay BDT histogram

// // Add a legend
// auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
// legend->AddEntry(DRN_ratio, "DRN Ratio", "l");
// legend->AddEntry(BDT_ratio, "BDT Ratio", "l");
// legend->Draw();

		// float Gen = matchedGenEnergy_->at(0);
		// float Reco_BDT = energy_ecal_mustache->at(0);

    //     recoDRNEnergy_->Fill(Reco_DRN);        matchedGenEnergy->Fill(Gen);        energy_ecal_mustache_->Fill(Reco_BDT);    

//     if (!recoDRNEnergy->empty()) {
//       float Reco_DRN = recoDRNEnergy->at(0);
  
//       recoDRNEnergy_->Fill(Reco_DRN);
//   } else {
//       cout << "Skipping event: Empty vector encountered!" << endl;
//   }

//   if (!matchedGenEnergy_->empty()) {
//     float Gen = matchedGenEnergy_->at(0);

//     matchedGenEnergy->Fill(Gen);
// } else {
//     cout << "Skipping event: Empty vector encountered!" << endl;
// }

// if (!energy_ecal_mustache->empty()) {
//   float Reco_BDT = energy_ecal_mustache->at(0);

//   energy_ecal_mustache_->Fill(Reco_BDT);
// } else {
//   cout << "Skipping event: Empty vector encountered!" << endl;
// }

	}

  DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
  BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
}

// DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
// BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
// ================================================= 2D Histograms =======================================================
// 	Al_eta_vs_Asl_eta->Fill(gen_a_lead_eta,gen_a_sublead_eta);
// 	Al_pt_vs_Asl_pt->Fill(gen_a_lead_pt,gen_a_sublead_pt);
// 	Al_p_vs_Asl_p->Fill(pA_lead_gen.P(),pA_sublead_gen.P());
// 	Al_mass_vs_Asl_mass ->Fill(gen_a_lead_mass,gen_a_sublead_mass);
// 	al_pho_eta1_eta2->Fill(gen_al_pho1_eta,gen_al_pho2_eta);
// 	asl_pho_eta1_eta2->Fill(gen_asl_pho1_eta,gen_asl_pho2_eta);
// 	al_pho_p1_p2 ->Fill(pPho1_al_gen.P(),pPho2_al_gen.P());
// 	asl_pho_p1_p2 ->Fill(pPho1_asl_gen.P(),pPho2_asl_gen.P());
// 	if(A_lead_EB){
// 	for (int i=0; i<A_flags->size();i++){
// 		if(A_flags->at(i)==0){
// 			float tot_unc_E=0;
// 		//cout<<"A_lead_eta"<<"\t"<<gen_a_lead_eta<<endl;
// 			for(int k=0; k< Pho_hit_x[i]->size();k++){
// 				a_eb_rechit_x->Fill(Pho_hit_x[i]->at(k));
// 				a_eb_rechit_y->Fill(Pho_hit_y[i]->at(k));
// 				a_eb_rechit_eta->Fill(Pho_hit_eta[i]->at(k));
// 		//cout<<Pho_hit_eta[i]->at(k)<<endl;

// 				a_eb_rechit_phi->Fill(Pho_hit_phi[i]->at(k));
// 				a_eb_hit_eta_phi->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k));
// 				a_eb_hit_xy->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k));
// 				a_eb_hit_eta_phi_En->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k),Pho_hit_E[i]->at(k));
// 				a_eb_hit_xy_En->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k),Pho_hit_E[i]->at(k));
// 			if(Pho_hit_frac[i]->at(k)<0){tot_unc_E = tot_unc_E + Pho_hit_E[i]->at(k);}
// 					}
// 			Tot_unc_E->Fill(tot_unc_E);
// 			Tot_unc_E_eb->Fill(tot_unc_E);
// 				}
// 			}
// 		}	
// 	if(A_sublead_EB){
// 	for (int i=0; i<A_flags->size();i++){
// 		if(A_flags->at(i)==1){
// 			float tot_unc_E = 0;
// 			for(int k=0; k< Pho_hit_x[i]->size();k++){
// 				a_eb_rechit_x->Fill(Pho_hit_x[i]->at(k));
// 				a_eb_rechit_y->Fill(Pho_hit_y[i]->at(k));
//                                 a_eb_rechit_eta->Fill(Pho_hit_eta[i]->at(k));
//                                 a_eb_rechit_phi->Fill(Pho_hit_phi[i]->at(k));
// 				a_eb_hit_eta_phi->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k));
//                                 a_eb_hit_xy->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k));
//                                 a_eb_hit_eta_phi_En->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k),Pho_hit_E[i]->at(k));
//                                 a_eb_hit_xy_En->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k),Pho_hit_E[i]->at(k));
// 			if(Pho_hit_frac[i]->at(k)<0){tot_unc_E = tot_unc_E + Pho_hit_E[i]->at(k);}
// 					}
// 			Tot_unc_E->Fill(tot_unc_E);
// 			Tot_unc_E_eb->Fill(tot_unc_E);
// 				}
// 			}
// 		}	
       
// 	if(A_lead_EE){
// 	for (int i=0; i<A_flags->size();i++){
// 		if(A_flags->at(i)==0){
// 			float tot_unc_E=0;
// //cout<<"A_lead_eta"<<"\t"<<gen_a_lead_eta<<endl;
// 			for(int k=0; k< Pho_hit_x[i]->size();k++){
// 				a_ee_rechit_x->Fill(Pho_hit_x[i]->at(k));
// 				a_ee_rechit_y->Fill(Pho_hit_y[i]->at(k));
// 				a_ee_rechit_eta->Fill(Pho_hit_eta[i]->at(k));
// //cout<<Pho_hit_eta[i]->at(k)<<endl;
// 				a_ee_rechit_phi->Fill(Pho_hit_phi[i]->at(k));
// 				a_ee_hit_eta_phi->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k));
// 				a_ee_hit_xy->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k));
// 				a_ee_hit_eta_phi_En->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k),Pho_hit_E[i]->at(k));
// 				a_ee_hit_xy_En->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k),Pho_hit_E[i]->at(k));
// 			if(Pho_hit_frac[i]->at(k)<0){tot_unc_E = tot_unc_E + Pho_hit_E[i]->at(k);}
// 					}
// 			Tot_unc_E->Fill(tot_unc_E);
// 			Tot_unc_E_ee->Fill(tot_unc_E);
// 				}
// 			}
// 		}	
// 	if(A_sublead_EE){
// 	for (int i=0; i<A_flags->size();i++){
// 		if(A_flags->at(i)==1){
// 		float tot_unc_E=0;
// 			for(int k=0; k< Pho_hit_x[i]->size();k++){
// 				a_ee_rechit_x->Fill(Pho_hit_x[i]->at(k));
// 				a_ee_rechit_y->Fill(Pho_hit_y[i]->at(k));
// 				a_ee_rechit_eta->Fill(Pho_hit_eta[i]->at(k));
// 				a_ee_rechit_phi->Fill(Pho_hit_phi[i]->at(k));
// 				a_ee_hit_xy->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k));
// 				a_ee_hit_eta_phi->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k));
// 				a_ee_hit_eta_phi_En->Fill(Pho_hit_eta[i]->at(k),Pho_hit_phi[i]->at(k),Pho_hit_E[i]->at(k));
//                                 a_ee_hit_xy_En->Fill(Pho_hit_x[i]->at(k),Pho_hit_y[i]->at(k),Pho_hit_E[i]->at(k));
// 		if(Pho_hit_frac[i]->at(k)<0){tot_unc_E = tot_unc_E + Pho_hit_E[i]->at(k);}
// 					}
// 			Tot_unc_E->Fill(tot_unc_E);
// 			Tot_unc_E_ee->Fill(tot_unc_E);
// 				}
// 			}
// 		}	
 
//  }
//  cout<<"Veto Counts:"<<"\t"<<veto_count<<endl;
// }
