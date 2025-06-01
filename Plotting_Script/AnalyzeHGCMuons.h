#ifndef ANALYZEHGCMuons_H
#define ANALYZEHGCMuons_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
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
class AnalyzeHGCMuons : public HGCNtupleVariables
{

public:
  AnalyzeHGCMuons(const TString &inputFileList = "foo.txt", const char *outFileName = "histo.root", const char *dataset = "data", const char *massP="0.1");
  ~AnalyzeHGCMuons();
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

  TH1F * DRN_ratio;
  TH1F * BDT_ratio;


 /* vector<TH2F*> TwoDHist;= {a_gen_mass_vs_pt,a_gen_eta_vs_phi, angle_vs_gamma_Ma_240_260,angle_vs_gamma_Ma_490_510,angle_vs_gamma_Ma_740_760, angle_vs_gamma_Ma_990_1010,
  deta_dphi_gen, EE_XY_occupancy, EE_XY_occu_En_weighed};*/
};
#endif

#ifdef ANALYZEHGCMuons_cxx

void AnalyzeHGCMuons::BookHistogram(const char *outFileName)
{

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();
  if (fillhist)
  { //Define histograms
   DRN_ratio      =   new       TH1F    ("DRN_corrbyGen"                ,     "E_DRN/E_Gen"       ,    100000   ,   0     ,      2.5);
   BDT_ratio     =   new       TH1F    ("BDT_corrbyGen"                ,     "E_BDT/E_Gen"       ,    100000   ,   0     ,      2.5);
    // energy_ecal_mustache_      =   new       TH1F    ("energy_ecal_mustache_"                ,     "BDT corrected energy"       ,    500   ,   0.5     ,      1.5);
    
// ECAL Hits

// 2D Histograms

  }
  }

AnalyzeHGCMuons::AnalyzeHGCMuons(const TString &inputFileList, const char *outFileName, const char *dataset, const char *massP)
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

  HGCNtupleVariables::Init(tree, tree2);

  BookHistogram(outFileName);
//cout<<massP<<endl;
  // treeout = new TTree("fordrn", "ForDRN");
  // treeout->Branch("Hit_X_Pho1", &M_HitX1);
  // treeout->Branch("EEFlag1", &M_EEFlag1);

}

Bool_t AnalyzeHGCMuons::FillChain(TChain *chain, TChain *chain2, const TString &inputFileList)
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

Long64_t AnalyzeHGCMuons::LoadTree(Long64_t entry)
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
void AnalyzeHGCMuons::DrawHist(TH1F *hist1,string name, string x_title="", string y_title="")
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
void AnalyzeHGCMuons::DrawHistLogy(TH1F *hist1,string name, string x_title="", string y_title="")
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
void AnalyzeHGCMuons::DrawHistTXT(TH2F *hist1, string name, string x_title="", string y_title="")
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
void AnalyzeHGCMuons::DrawHist(TH2F *hist1, string name, string x_title="", string y_title="")
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
void AnalyzeHGCMuons::DrawHistNoStat(TH2F *hist1, string name, string x_title="", string y_title="")
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


void AnalyzeHGCMuons::graphSide(TH1F *hist1, TH1F *hist2, string name,string x_title="", string y_title="")
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
void AnalyzeHGCMuons::graphSide(TH2F *hist1, TH2F *hist2, string name,string x_title, string y_title)
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

void AnalyzeHGCMuons::graphOverlay(TH1F *hist1, TH1F *hist2, TH1F *hist4, TH1F *hist5, string name,string x_title="",string y_title="")
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
void AnalyzeHGCMuons::graphOverlayLogy(TH1F *hist1, TH1F *hist2, TH1F *hist4, TH1F *hist5, string name,string x_title="",string y_title="")
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
void AnalyzeHGCMuons::graphOverlay(TH1F *hist1, TH1F *hist2, string name,string x_title ="",string y_title="")
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

void AnalyzeHGCMuons::graphOverlay(TH1F *hist1, TH1F *hist2, TH1F *hist4,  string name,string x_title="",string y_title="")
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


void AnalyzeHGCMuons::graphHist(TH2F *hist, char title[], char name[])
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

AnalyzeHGCMuons::~AnalyzeHGCMuons()
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
