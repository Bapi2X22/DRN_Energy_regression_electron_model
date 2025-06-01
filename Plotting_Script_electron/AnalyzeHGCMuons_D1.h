#ifndef ANALYZEHGCMuons_D1_H
#define ANALYZEHGCMuons_D1_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables_D1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TText.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFrame.h"
#include "TFile.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TVirtualHistPainter.h"
#include "TPaveStats.h"
#include "TLatex.h"

void PaintOverflow(TH1 *h)
{
  // This function paint the histogram h with an extra bin for overflows

  const char *name = h->GetName();
  const char *title = h->GetTitle();
  Int_t nx = h->GetNbinsX() + 1;
  Double_t x1 = h->GetBinLowEdge(1);
  Double_t bw = h->GetBinWidth(nx);
  Double_t x2 = h->GetBinLowEdge(nx) + bw;

  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F("rebinned", title, nx, x1, x2);

  // Fill the new hitogram including the extra bin for overflows
  for (Int_t i = 1; i <= nx; i++)
  {
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  }

  // Fill the underflows
  htmp->Fill(x1 - 1, h->GetBinContent(0));

  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());

  // Draw the temporary histogram
  htmp->Draw();
  TText *t = new TText(x2 - bw / 2, h->GetBinContent(nx), "Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.03);
  ;
  t->Draw();
}
class AnalyzeHGCMuons_D1 : public HGCNtupleVariables_D1
{

public:
  AnalyzeHGCMuons_D1(const TString &inputFileList = "foo.txt", const char *outFileName = "histo.root", const char *dataset = "data", const char *massP="0.1");
  ~AnalyzeHGCMuons_D1();
  Bool_t FillChain(TChain *chain, TChain *chain2, const TString &inputFileList);
  void graphHist(TH2F *hist, char title[], char name[]);
  Long64_t LoadTree(Long64_t entry);
  void EventLoop(const char *);
  void BookHistogram(const char *);
  void graphSide(TH1F *hist1, TH1F *hist2, string name,string x_title, string y_title);
  void graphSide(TH2F *hist1, TH2F *hist2, string name, string x_title, string y_title);
  void graphOverlay(TH2F *hist1, TH2F *hist2, string name,string x_title,string y_title);
  void graphOverlay(TH1F *hist1, TH1F *hist2, string name,string x_title,string y_title);
  void graphOverlay(TH1F *hist1, TH1F *hist2, TH1F *hist4,string name,string x_title,string y_title);
  void graphOverlay(TH1F *hist1, TH1F *hist2,TH1F *hist4,TH1F *hist5, string name,string x_title,string y_title);
  void graphOverlayLogy(TH1F *hist1, TH1F *hist2,TH1F *hist4,TH1F *hist5, string name,string x_title,string y_title);
  void DrawHist(TH1F *hist1, string name, string x_title, string y_title);
  void DrawHistLogy(TH1F *hist1, string name, string x_title, string y_title);
  void DrawHist(TH2F *hist1, string name, string x_title, string y_title);
  void DrawHistTXT(TH2F *hits, string name, string x_title, string y_title);
  void DrawHistNoStat(TH2F *hist1, string name, string x_title, string y_title);
  int fillhist = 1;
  /*int ctr_twogen = 0;
  int ctr_mcut = 0;
  int ctr_tworeco = 0;
  int ctr_onereco = 0;
  int ctr_onedr = 0;
  int ctr_twodr = 0;*/
//  string mass = "1.8";
  string gamma = "";
  string mass = "1.8";

  TFile *oFile;
  TDirectory *d_Layer1;

  // TTree *treeout;
  // double M_m_gen;
  // vector<float> M_pt_reco;

  // vector<float> *M_HitX1;
  // bool M_EEFlag1;
  // bool M_EEFlag2;
  // vector<float> *M_HitY1;
  // vector<float> *M_HitZ1;
  // vector<float> *M_HitEta1;
  // vector<float> *M_HitPhi1;


  // double M_EBM;
  // double M_angle_gen;
  // int M_event;
  // int M_lumi;
  // int M_run;
  //  float m_gen;
  //Declare Historgrams
  //1D histograms
  // TH1F * recoDRNEnergy_;
  // TH1F * matchedGenEnergy;
  // TH1F * energy_ecal_mustache_;
  int drnCloserCount = 0;
  int bdtCloserCount = 0;
  TH1F * DRN_ratio;
  TH1F * BDT_ratio;
  TH1F * mee_BDT;
  TH1F * mee_DRN;
  TH1F * Num_rechits;
  TH1F * Energy_rechits;
  TH1F * Sum_En_rechits;
  TH1F * R9;
  TH1F * R9_EB;
  TH1F * R9_EE;
  TH1F * SigIEIE;
  TH1F * SigIEIE_EB;
  TH1F * SigIEIE_EE;
  TH1F * HitNoise;
  std::vector<TH1D*> pT_histos_BDT;
  std::vector<TH1D*> pT_histos_DRN;
  std::vector<TH1D*> R9_histos_BDT;
  std::vector<TH1D*> R9_histos_DRN;
  std::vector<TH1D*> eta_histos_BDT_EE;
  std::vector<TH1D*> eta_histos_DRN_EE;
  std::vector<TH1D*> eta_histos_BDT_EB;
  std::vector<TH1D*> eta_histos_DRN_EB;
  // std::vector<std::pair<double, double>> pT_bins = {
  //       {5, 10},
  //       {10, 50},
  //       {50, 200},
  //       {200, 500},
  //       {500, 1000},
  //       {1000, 1500},
  //       {1500, 5000}
  //   };
 // Define pT bin ranges from 0 to 600 in steps of 50
// std::vector<std::pair<double, double>> pT_bins = {
//     {0, 50}, {50, 100}, {100, 150}, {150, 200},
//     {200, 250}, {250, 300}, {300, 350}, {350, 400},
//     {400, 450}, {450, 500}, {500, 550}, {550, 600}
// };

std::vector<std::pair<int, int>> pT_bins = {
    {0, 20}, {20, 40}, {40, 60}, {60, 80},
    {80, 100}, {100, 120}, {120, 140}, {140, 160},
    {160, 180}, {180, 200}, {200, 220}, {220, 240},
    {240, 260}, {260, 280}, {280, 300}, {300, 320},
    {320, 340}, {340, 360}, {360, 380}, {380, 400},
    {400, 420}, {420, 440}, {440, 460}, {460, 480},
    {480, 500}, {500, 520}, {520, 540}, {540, 560},
    {560, 580}, {580, 600}
};


// std::vector<std::pair<double, double>> R9_bins = {
//     {0, 0.75}, {0.75, 0.82}, {0.82, 0.85}, {0.85, 0.87},
// {0.87, 0.89}, {0.89, 0.9}, {0.91, 0.92}, {0.92, 0.93},
//     {0.93, 0.94}, {0.94, 0.95}, {0.95, 0.96}, {0.96, 0.97}, {0.97, 0.98},{0.98, 0.99},{0.99, 1.0}};

std::vector<std::pair<double, double>> R9_bins = {
    {0, 0.7}, {0.7, 0.8}, {0.8, 0.82}, {0.82, 0.84},
{0.84, 0.86}, {0.86, 0.88}, {0.88, 0.9}, {0.9, 0.91},
    {0.91, 0.92}, {0.92, 0.93}, {0.93, 0.94}, {0.94, 0.95}, {0.95, 0.96},{0.96, 0.97},{0.97, 0.98}, {0.98, 0.99}, {0.99, 1.0}, {1.0, 1.5}};

std::vector<std::pair<double, double>> EE_eta_bins = {
    {1.566, 1.7}, {1.7, 1.8}, {1.8, 1.9}, {1.9, 2.0},
{2.0, 2.1}, {2.1, 2.2}, {2.2, 2.3}, {2.3, 2.4},
    {2.4, 2.5}};

std::vector<std::pair<double, double>> EB_eta_bins = {
    {0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4},
{0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8},
    {0.8, 0.9}, {0.9, 1.0}, {1.0, 1.1}, {1.1, 1.2}, {1.2, 1.3},{1.3, 1.442}};



 /* vector<TH2F*> TwoDHist;= {a_gen_mass_vs_pt,a_gen_eta_vs_phi, angle_vs_gamma_Ma_240_260,angle_vs_gamma_Ma_490_510,angle_vs_gamma_Ma_740_760, angle_vs_gamma_Ma_990_1010,
  deta_dphi_gen, EE_XY_occupancy, EE_XY_occu_En_weighed};*/
};
#endif

#ifdef ANALYZEHGCMuons_D1_cxx

void AnalyzeHGCMuons_D1::BookHistogram(const char *outFileName)
{

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();
  if (fillhist)
  { //Define histograms
   DRN_ratio      =   new       TH1F    ("DRN_corrbyGen"                ,     "E_DRN/E_Gen"       ,    100000   ,   0     ,      2);
   DRN_ratio->Sumw2();
   BDT_ratio     =   new       TH1F    ("BDT_corrbyGen"                ,     "E_BDT/E_Gen"       ,    100000   ,   0     ,      2);
   BDT_ratio->Sumw2();
   mee_BDT = new TH1F    ("m_ee_BDT"                ,     "invarient_mass_BDT"       ,    100000   ,   0     ,      200);
   mee_BDT->Sumw2();
   mee_DRN = new TH1F    ("m_ee_DRN"                ,     "invarient_mass_DRN"       ,    100000   ,   0     ,      200);
   mee_DRN->Sumw2();
   Num_rechits = new       TH1F    ("num_rechits"                ,     "Number_of_rechits"       ,    1000   ,   0     ,      200);
   Num_rechits->Sumw2();
   Energy_rechits = new      TH1F    ("Energy_rechits"                ,     "Energy_of_rechits"       ,    10000   ,   0     ,      100);
   Energy_rechits->Sumw2();
   Sum_En_rechits = new      TH1F    ("sum_En_rechits"                ,     "Sum _of_Energy_of_rechits"       ,    10000   ,   0     ,      1000);
   Sum_En_rechits->Sumw2();
   R9 = new TH1F ("R9", "R9", 1000., 0, 2);
   R9->Sumw2();
   R9_EB = new TH1F ("R9_EB", "R9_EB", 1000, 0, 2);
   R9_EB->Sumw2();
   R9_EE = new TH1F ("R9_EE", "R9_EE", 1000, 0, 2);
   R9_EE->Sumw2();
   SigIEIE = new TH1F ("SigIEIE", "SigIEIE", 1000, 0, 0.2);
   SigIEIE->Sumw2();
   SigIEIE_EB = new TH1F ("SigIEIE_EB", "SigIEIE_EB", 1000, 0, 0.2);
   SigIEIE_EB->Sumw2();
   SigIEIE_EE = new TH1F ("SigIEIE_EE", "SigIEIE_EE", 1000, 0, 0.2);
   SigIEIE_EE->Sumw2();
   HitNoise = new TH1F ("HitNoise", "", 10000,0,10);
   HitNoise->Sumw2();

// Create histograms for each pT range with uniform bin width
   // for (const auto& bin : pT_bins) {
   //      std::string hist_name_B = "pT_bin_B_" + std::to_string(int(bin.first)) + "_" + std::to_string(int(bin.second));
   //      std::string hist_name_D = "pT_bin_D_" + std::to_string(int(bin.first)) + "_" + std::to_string(int(bin.second));
        
   //      // Define the number of bins uniformly within the range
   //      int num_bins = 100000;  // Change this value if you want a different bin granularity
   //      pT_histos_BDT.push_back(new TH1D(hist_name_B.c_str(), hist_name_B.c_str(), num_bins, 0, 2));
   //      pT_histos_DRN.push_back(new TH1D(hist_name_D.c_str(), hist_name_D.c_str(), num_bins, 0, 2));
   //  };

   //  for (const auto& bin : R9_bins) {
   //      std::string hist_name_B = "R9_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
   //      std::string hist_name_D = "R9_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
        
   //      // Define the number of bins uniformly within the range
   //      int num_bins = 100000;  // Change this value if you want a different bin granularity
   //      R9_histos_BDT.push_back(new TH1D(hist_name_B.c_str(), hist_name_B.c_str(), num_bins, 0, 2));
   //      R9_histos_DRN.push_back(new TH1D(hist_name_D.c_str(), hist_name_D.c_str(), num_bins, 0, 2));
   //  };
for (const auto& bin : pT_bins) {
    std::string hist_name_B = "pT_bin_B_" + std::to_string(int(bin.first)) + "_" + std::to_string(int(bin.second));
    std::string hist_name_D = "pT_bin_D_" + std::to_string(int(bin.first)) + "_" + std::to_string(int(bin.second));

    int num_bins = 100000;

    TH1D* hB = new TH1D(hist_name_B.c_str(), hist_name_B.c_str(), num_bins, 0, 2);
    hB->Sumw2();
    pT_histos_BDT.push_back(hB);

    TH1D* hD = new TH1D(hist_name_D.c_str(), hist_name_D.c_str(), num_bins, 0, 2);
    hD->Sumw2();
    pT_histos_DRN.push_back(hD);
}

for (const auto& bin : R9_bins) {
    std::string hist_name_B = "R9_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
    std::string hist_name_D = "R9_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);

    int num_bins = 100000;

    TH1D* hB = new TH1D(hist_name_B.c_str(), hist_name_B.c_str(), num_bins, 0, 2);
    hB->Sumw2();
    R9_histos_BDT.push_back(hB);

    TH1D* hD = new TH1D(hist_name_D.c_str(), hist_name_D.c_str(), num_bins, 0, 2);
    hD->Sumw2();
    R9_histos_DRN.push_back(hD);
}

for (const auto& bin : EE_eta_bins) {
    std::string hist_name_B = "eta_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
    std::string hist_name_D = "eta_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);

    int num_bins = 100000;

    TH1D* hB = new TH1D(hist_name_B.c_str(), hist_name_B.c_str(), num_bins, 0, 2);
    hB->Sumw2();
    eta_histos_BDT_EE.push_back(hB);

    TH1D* hD = new TH1D(hist_name_D.c_str(), hist_name_D.c_str(), num_bins, 0, 2);
    hD->Sumw2();
    eta_histos_DRN_EE.push_back(hD);
}

for (const auto& bin : EB_eta_bins) {
    std::string hist_name_B = "eta_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
    std::string hist_name_D = "eta_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);

    int num_bins = 100000;

    TH1D* hB = new TH1D(hist_name_B.c_str(), hist_name_B.c_str(), num_bins, 0, 2);
    hB->Sumw2();
    eta_histos_BDT_EB.push_back(hB);

    TH1D* hD = new TH1D(hist_name_D.c_str(), hist_name_D.c_str(), num_bins, 0, 2);
    hD->Sumw2();
    eta_histos_DRN_EB.push_back(hD);
}

 
    // energy_ecal_mustache_      =   new       TH1F    ("energy_ecal_mustache_"                ,     "BDT corrected energy"       ,    500   ,   0.5     ,      1.5);
    
// ECAL Hits

// 2D Histograms

  }
  }

AnalyzeHGCMuons_D1::AnalyzeHGCMuons_D1(const TString &inputFileList, const char *outFileName, const char *dataset, const char *massP)
{
    TChain *tree = new TChain("nTuplelize/T");
  TChain *tree2;
  mass = string(massP);

  if (!FillChain(tree, tree2, inputFileList))
    {std::cerr << "Cannot get the tree " << std::endl;
  }
  else
  {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;  }

  HGCNtupleVariables_D1::Init(tree, tree2);

  BookHistogram(outFileName);
//cout<<massP<<endl;
  // treeout = new TTree("fordrn", "ForDRN");
  // treeout->Branch("Hit_X_Pho1", &M_HitX1);
  // treeout->Branch("EEFlag1", &M_EEFlag1);

}

Bool_t AnalyzeHGCMuons_D1::FillChain(TChain *chain, TChain *chain2, const TString &inputFileList)
{

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if (!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while (1)
  {
    infile >> buffer;
    if (!infile.good())
      break;
    // std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
    // chain2->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeHGCMuons_D1::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (!fChain->InheritsFrom(TChain::Class()))
    return centry;
  TChain *chain = (TChain *)fChain;
  if (chain->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  // Modified By Chirayu
  return centry;
  // End Modified

  if (!fChain2)
    return -5;
  Long64_t centry2 = fChain2->LoadTree(entry);
  if (centry2 < 0)
    return centry2;
  if (!fChain2->InheritsFrom(TChain::Class()))
    return centry2;
  TChain *chain2 = (TChain *)fChain2;
  if (chain2->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (centry == centry2)
    return centry;
  else
    return -1;
}
void AnalyzeHGCMuons_D1::DrawHist(TH1F *hist1,string name, string x_title="", string y_title="")
{
string Name = name + ";" + x_title + ";" + y_title;
TCanvas *C= new TCanvas("C","C",1000,1000);
C->cd();
hist1->Draw("hist");
hist1->SetTitle(Name.c_str());
hist1->GetXaxis()->CenterTitle(true);
hist1->GetYaxis()->CenterTitle(true);
hist1->Write();
C->Update();
TPaveStats *st1 = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");	
st1->SetX1NDC(.77); st1->SetX2NDC(.9);st1->SetY1NDC(.9); st1->SetY2NDC(.77);
st1->Draw();
C->Modified();
C->Write();
C->SaveAs((name + string(".png")).c_str());
C->SaveAs((name + string(".pdf")).c_str());

}
void AnalyzeHGCMuons_D1::DrawHistLogy(TH1F *hist1,string name, string x_title="", string y_title="")
{
string Name = name + ";" + x_title + ";" + y_title;
TCanvas *C= new TCanvas("C","C",1000,1000);
C->cd();
hist1->Draw("hist");
hist1->SetTitle(Name.c_str());
hist1->GetXaxis()->CenterTitle(true);
hist1->GetYaxis()->CenterTitle(true);
hist1->Write();
C->SetLogy();
C->Update();
TPaveStats *st1 = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");	
st1->SetX1NDC(.77); st1->SetX2NDC(.9);st1->SetY1NDC(.9); st1->SetY2NDC(.77);
st1->Draw();
C->Modified();
C->Write();
C->SaveAs((name + string(".png")).c_str());
C->SaveAs((name + string(".pdf")).c_str());

}
void AnalyzeHGCMuons_D1::DrawHistTXT(TH2F *hist1, string name, string x_title="", string y_title="")
{
string Name = name + ";" + x_title + ";" + y_title;
TCanvas *C1= new TCanvas (name.c_str(),"",1200,1000);
C1->SetRightMargin(0.15);
hist1->Draw("COLZ+TEXT");
hist1->SetTitle(Name.c_str());
hist1->GetXaxis()->CenterTitle(true);
hist1->GetYaxis()->CenterTitle(true);
C1->Update();
hist1->Write();
hist1->SetStats(0);
C1->Modified();
C1->Write();
C1->SaveAs((name + string(".png")).c_str());
C1->SaveAs((name + string(".pdf")).c_str());

}
void AnalyzeHGCMuons_D1::DrawHist(TH2F *hist1, string name, string x_title="", string y_title="")
{
string Name = name + ";" + x_title + ";" + y_title;
TCanvas *C1= new TCanvas (name.c_str(),"",1200,1000);
C1->SetRightMargin(0.15);
hist1->Draw("COLZ");
hist1->SetTitle(Name.c_str());
hist1->GetXaxis()->CenterTitle(true);
hist1->GetYaxis()->CenterTitle(true);
hist1->Write();
C1->Modified();

C1->Update();
C1->Write();
C1->SaveAs((name + string(".png")).c_str());
C1->SaveAs((name + string(".pdf")).c_str());

}
void AnalyzeHGCMuons_D1::DrawHistNoStat(TH2F *hist1, string name, string x_title="", string y_title="")
{
string Name = name + ";" + x_title + ";" + y_title;
TCanvas *C1= new TCanvas (name.c_str(),"",1200,1000);
C1->SetRightMargin(0.15);
hist1->Draw("COLZ");
hist1->SetTitle(Name.c_str());
hist1->GetXaxis()->CenterTitle(true);
hist1->GetYaxis()->CenterTitle(true);
C1->Update();
hist1->SetStats(0);
hist1->Write();
C1->Modified();

C1->Write();
C1->SaveAs((name + string(".png")).c_str());
C1->SaveAs((name + string(".pdf")).c_str());

}


void AnalyzeHGCMuons_D1::graphSide(TH1F *hist1, TH1F *hist2, string name,string x_title="", string y_title="")
{
  TCanvas *c1 = new TCanvas(name.c_str(), "", 200, 10, 700, 500);
  c1->Divide(2, 1);
  c1->cd(1);
  hist1->Draw("colz");
  hist1->GetXaxis()->SetTitle(x_title.c_str());
  hist1->GetYaxis()->SetTitle(y_title.c_str());
  hist1->GetXaxis()->CenterTitle(true);
  hist1->GetYaxis()->CenterTitle(true);
  c1->cd(2);
  hist2->Draw("colz");
  hist2->GetXaxis()->SetTitle(x_title.c_str());
  hist2->GetYaxis()->SetTitle(y_title.c_str());
  hist2->GetXaxis()->CenterTitle(true);
  hist2->GetYaxis()->CenterTitle(true);
  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
  c1->Write();
  c1->SaveAs((name + string(".png")).c_str());
 c1->SaveAs((name + string(".pdf")).c_str());
}
void AnalyzeHGCMuons_D1::graphSide(TH2F *hist1, TH2F *hist2, string name,string x_title, string y_title)
{
  TCanvas *c1 = new TCanvas(name.c_str(), "", 200, 10, 700, 500);
  c1->Divide(2, 1);
  c1->cd(1);
  hist1->Draw("colz");
  hist1->GetXaxis()->SetTitle(x_title.c_str());
  hist1->GetYaxis()->SetTitle(x_title.c_str());
  hist1->GetXaxis()->CenterTitle(true);
  hist1->GetYaxis()->CenterTitle(true);

  c1->cd(2);
  hist2->Draw("colz");
  hist1->GetXaxis()->SetTitle(x_title.c_str());
  hist1->GetYaxis()->SetTitle(x_title.c_str());
  hist1->GetXaxis()->CenterTitle(true);
  hist1->GetYaxis()->CenterTitle(true);

  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
  c1->Write();
  c1->SaveAs((name + string(".png")).c_str());
  c1->SaveAs((name + string(".png")).c_str());

}

void AnalyzeHGCMuons_D1::graphOverlay(TH1F *hist1, TH1F *hist2, TH1F *hist4, TH1F *hist5, string name,string x_title="",string y_title="")
{ 
  string Name=name + ";"+ x_title + ";" + y_title;
  gStyle->SetOptStat(1111);
  THStack *hist3 = new THStack(name.c_str(), Name.c_str());
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  hist4->SetLineColor(kGreen);
  hist5->SetLineColor(kBlack);
  hist3->Add(hist1,"sames");
  hist3->Add(hist2,"sames");
  hist3->Add(hist4,"sames");
  hist3->Add(hist5,"sames");
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  c3->cd();
  hist3->Draw("NOSTACK");
    hist3->Write();
  c3->Update();
  TPaveStats *st1 = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st2 = (TPaveStats*)hist2->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st4 = (TPaveStats*)hist4->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st5 = (TPaveStats*)hist5->GetListOfFunctions()->FindObject("stats");	
st1->SetX1NDC(.7); st1->SetX2NDC(.9);st1->SetY1NDC(.7); st1->SetY2NDC(.6);
st2->SetX1NDC(.7); st2->SetX2NDC(.9);st2->SetY1NDC(.6); st2->SetY2NDC(.5);
st4->SetX1NDC(.7); st4->SetX2NDC(.9);st4->SetY1NDC(.5); st4->SetY2NDC(.4);
st5->SetX1NDC(.7); st5->SetX2NDC(.9);st5->SetY1NDC(.4); st5->SetY2NDC(.3);
st1->Draw();st2->Draw();
st4->Draw();st5->Draw();
  c3->Modified();
  c3->BuildLegend(0.7, 0.7, 0.9, 0.9, "");
  
  c3->Write();
 
  c3->SaveAs((name + string(".pdf")).c_str());
  c3->SaveAs((name + string(".png")).c_str());
}
void AnalyzeHGCMuons_D1::graphOverlayLogy(TH1F *hist1, TH1F *hist2, TH1F *hist4, TH1F *hist5, string name,string x_title="",string y_title="")
{ 
  string Name=name + ";"+ x_title + ";" + y_title;
  gStyle->SetOptStat(1111);
  THStack *hist3 = new THStack(name.c_str(), Name.c_str());
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  hist4->SetLineColor(kGreen);
  hist5->SetLineColor(kBlack);
  hist3->Add(hist1,"sames");
  hist3->Add(hist2,"sames");
  hist3->Add(hist4,"sames");
  hist3->Add(hist5,"sames");
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  c3->cd();
  c3->SetLogy();
  hist3->Draw("NOSTACK");
    hist3->Write();
  c3->Update();
  TPaveStats *st1 = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st2 = (TPaveStats*)hist2->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st4 = (TPaveStats*)hist4->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st5 = (TPaveStats*)hist5->GetListOfFunctions()->FindObject("stats");	
st1->SetX1NDC(.7); st1->SetX2NDC(.9);st1->SetY1NDC(.7); st1->SetY2NDC(.6);
st2->SetX1NDC(.7); st2->SetX2NDC(.9);st2->SetY1NDC(.6); st2->SetY2NDC(.5);
st4->SetX1NDC(.7); st4->SetX2NDC(.9);st4->SetY1NDC(.5); st4->SetY2NDC(.4);
st5->SetX1NDC(.7); st5->SetX2NDC(.9);st5->SetY1NDC(.4); st5->SetY2NDC(.3);
st1->Draw();st2->Draw();
st4->Draw();st5->Draw();
  c3->Modified();
  c3->BuildLegend(0.7, 0.7, 0.9, 0.9, "");
  
  c3->Write();
 
  c3->SaveAs((name + string(".pdf")).c_str());
  c3->SaveAs((name + string(".png")).c_str());
}
void AnalyzeHGCMuons_D1::graphOverlay(TH1F *hist1, TH1F *hist2, string name,string x_title ="",string y_title="")
{ 
  string Name =name + ";"+ x_title + ";" + y_title;
  gStyle->SetOptStat(1111);
  THStack *hist3 = new THStack(name.c_str(), name.c_str());
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  hist3->Add(hist1,"sames");
  hist3->Add(hist2,"sames");
   TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
 // c3->cd();
  hist3->Draw("NOSTACK");
  hist3->GetXaxis()->SetTitle(x_title.c_str());
  hist3->GetYaxis()->SetTitle(y_title.c_str());
   hist3->Write();
  c3->Update();
  TPaveStats *st1 = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st2 = (TPaveStats*)hist2->GetListOfFunctions()->FindObject("stats");	
st1->SetX1NDC(.7); st1->SetX2NDC(.9);st1->SetY1NDC(.7); st1->SetY2NDC(.6);
st2->SetX1NDC(.7); st2->SetX2NDC(.9);st2->SetY1NDC(.6); st2->SetY2NDC(.5);
st1->Draw();st2->Draw();
  c3->Modified();
  c3->BuildLegend(0.7, 0.7, 0.9, 0.9, "");
 
  c3->Write();
 
  c3->SaveAs((name + string(".pdf")).c_str());
  c3->SaveAs((name + string(".png")).c_str());
}

void AnalyzeHGCMuons_D1::graphOverlay(TH1F *hist1, TH1F *hist2, TH1F *hist4,  string name,string x_title="",string y_title="")
{
  string Name=name + ";"+ x_title + ";" + y_title;
   gStyle->SetOptStat(1111);
  //c3->cd();
  THStack *hist3 = new THStack(name.c_str(), Name.c_str());
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  hist4->SetLineColor(kGreen);
  
  hist3->Add(hist1,"sames");
  hist3->Add(hist2,"sames");
  hist3->Add(hist4,"sames");
  
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  hist3->Draw("NOSTACK");
    hist3->Write();
  c3->Update();
  TPaveStats *st1 = (TPaveStats*)hist1->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st2 = (TPaveStats*)hist2->GetListOfFunctions()->FindObject("stats");	
  TPaveStats *st4 = (TPaveStats*)hist4->GetListOfFunctions()->FindObject("stats");	
  
st1->SetX1NDC(.7); st1->SetX2NDC(.9);st1->SetY1NDC(.7); st1->SetY2NDC(.6);
st2->SetX1NDC(.7); st2->SetX2NDC(.9);st2->SetY1NDC(.6); st1->SetY2NDC(.5);
st4->SetX1NDC(.7); st4->SetX2NDC(.9);st4->SetY1NDC(.5); st1->SetY2NDC(.4);

st1->Draw();st2->Draw();
st4->Draw();
  c3->Modified();
  c3->BuildLegend(0.7, 0.7, 0.9, 0.9, "");
  
  c3->Write();
 
  c3->SaveAs((name + string(".pdf")).c_str());
  c3->SaveAs((name + string(".png")).c_str());
}


void AnalyzeHGCMuons_D1::graphHist(TH2F *hist, char title[], char name[])
{

  TCanvas *c1 = new TCanvas(name, title, 200, 10, 700, 500);
  TGraph *gr1 = new TGraph();
  for (int i = 0; i < hist->GetNbinsX(); i++)
  {
    for (int j = 0; j < hist->GetNbinsY(); j++)
    {
      double x = hist->GetXaxis()->GetBinCenter(i);
      double y = hist->GetYaxis()->GetBinCenter(j);
      gr1->SetPoint(gr1->GetN(), x, y);
    }
  }
  gr1->SetLineColor(2);
  gr1->SetLineWidth(4);
  gr1->SetMarkerStyle(0);
  gr1->SetTitle(title);
  gr1->GetXaxis()->SetTitle("Eta");
  gr1->GetYaxis()->SetTitle("Pt");
  gr1->Draw("ACP");
  gr1->Write(name);

  // TCanvas::Update() draws the frame, after which one can change it
  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
}
void DivideHistogramZValues(TH2F *hist, Double_t constant)
{
  Int_t nBinsX = hist->GetNbinsX();
  Int_t nBinsY = hist->GetNbinsY();

  // Iterate over each bin and divide the z-value by the constant
  for (Int_t i = 1; i <= nBinsX; i++)
  {
    for (Int_t j = 1; j <= nBinsY; j++)
    {
      Double_t binContent = hist->GetBinContent(i, j);
      Double_t newContent = binContent / constant;
      hist->SetBinContent(i, j, newContent);
    }
  }
}

AnalyzeHGCMuons_D1::~AnalyzeHGCMuons_D1()
{

  fillhist = 1;
  if (fillhist)
  {
 
 }


  // if (!fChain || !fChain2) return;
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
  // delete fChain2->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  // fChain->CloneTree(-1,"");
  oFile->Close();
}

#endif
// TODO: Pt condition in dR matching (later), Plots vs mass using jupyter notebook, Publish , Write explainations( tomorrow)
