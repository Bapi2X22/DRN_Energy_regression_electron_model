int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};                                                                   
int marker_style[12] = {8,22,34,47,33,22,34,39,41,43,45,23};
int line_color[9] = {kBlue,kRed,kGreen+2,kYellow+1,kViolet+2,kCyan+2,kGray+2,kMagenta,kBlue+2};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
vector<int> col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
// void decorate(TH1F*,int);
// void decorate(TH1F* hist,int i){
//   hist->SetLineWidth(3);                                                                                                                        
// }
// void setLastBinAsOverFlow(TH1F*);
// TH1F* setMyRange(TH1F*,double,double);

// TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
//   //call it after setting last bin as overflow                                                                                                                               
//   double err=0;
//   if(xHigh > 13000) return h1;
//   if(xLow < -13000) return h1;                                                                                                                     
//   int nMax=h1->FindBin(xHigh);
//   h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
//   h1->SetBinError(nMax,err);                                                                                                                 
//   for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
//     h1->SetBinContent(i,0);
//     h1->SetBinError(i,0);
                                                                                                                         
//   }
//   h1->GetXaxis()->SetRangeUser(xLow,xHigh);
//   cout<<xLow<<"\t"<<xHigh<<"\t"<<"set range"<<endl;
//    return h1;
// }

// TH1F* DrawOverflow(TH1F*);
// TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){

//    UInt_t nx    = h->GetNbinsX()+1;
//    Double_t *xbins= new Double_t[nx+1];
//    for (UInt_t i=0;i<nx;i++)
//      xbins[i]=h->GetBinLowEdge(i+1);
//    xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
//    char *tempName= new char[strlen(h->GetName())+10];
//    sprintf(tempName,"%swtOverFlow",h->GetName());
//    h->GetXaxis()->SetLimits(xmin,xrange);
//    // Book a temporary histogram having ab extra bin for overflows
//    TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
//    htmp->GetXaxis()->SetRange(xmin,xrange);
//    // Reset the axis labels
//    htmp->SetXTitle(h->GetXaxis()->GetTitle());
//    htmp->SetYTitle(h->GetYaxis()->GetTitle());
//    // Fill the new hitogram including the extra bin for overflows
//    for (UInt_t i=1; i<=nx; i++)
//      htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
//    // Fill the underflows
//    htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   
//    htmp->SetEntries(h->GetEntries());
   
//    return htmp;
// }
// void setLastBinAsOverFlow(TH1F* h_hist){
//   double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
//   double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

//   if(lastBinCt!=0 && overflCt!=0)
//     lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

//   else if(lastBinCt==0 && overflCt!=0)
//     lastBinErr = overflErr;
//   else if(lastBinCt!=0 && overflCt==0)
//     lastBinErr = lastBinErr;
//   else lastBinErr=0;

//   lastBinCt = lastBinCt+overflCt;
//   h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
//   h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
//   cout<<lastBinCt<<"\t"<<"Last bin values"<<endl;

// }
#include "TF1.h"
#include "TMath.h"
#include <cmath>

// Define the Cruijff function
double Cruijff(double *x, double *par) {
    double dx = x[0] - par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double sigma = (sigmaL + sigmaR) / 2;
    double alphaL = par[3];
    double alphaR = par[4];
    
    if (dx < 0)
        return par[5] * TMath::Exp(-dx * dx / (2 * sigmaL * sigmaL + alphaL * dx * dx));
    else
        return par[5] * TMath::Exp(-dx * dx / (2 * sigmaR * sigmaR + alphaR * dx * dx));
}

const int nfiles=100;                                                                                                                                                             
TFile *f[nfiles];

void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", char const *xlabel="", char const *ylabel="", int rebin=-1, double ymin=0, double ymax=0, double xmin=-1, double xmax=-1, char const *leg_head="", bool normalize=false, bool log_flag=false, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"})
{   normalize=true;
    TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name, 950, 850);
    canvas_n1->SetLeftMargin(0.135);
    canvas_n1->SetRightMargin(0.035);
    canvas_n1->SetTopMargin(0.04);
    canvas_n1->SetBottomMargin(0.1);
    gStyle->SetOptStat(0);

    TLegend *legend = new TLegend(0.21, 0.82, 0.65, 0.95);
    legend->SetTextSize(0.035);
    legend->SetLineColor(kWhite);
    legend->SetNColumns(2);
    legend->SetHeader(title);
    double MaxY=0;
    for(int i = 0; i < (int)hist.size(); i++) {
        hist.at(i)->Rebin(rebin);
        hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);
        if(hist.at(i)->GetMaximum()/hist.at(i)->Integral() > MaxY){MaxY=hist.at(i)->GetMaximum()/hist.at(i)->Integral();}
                                              }
 
 
    vector<vector<TLatex>> latexVec(hist.size(),vector<TLatex>(3));
    vector<TF1*> Fit_func(hist.size());
    for(int i = 0; i < (int)hist.size(); i++) {
        if(normalize && hist.at(i)->Integral() > 0) {
            hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
            hist.at(i)->GetYaxis()->SetTitle("Normalized");
        } else {
            hist.at(i)->GetYaxis()->SetTitle("Entries");
        }

        hist.at(i)->SetLineWidth(2);
        //hist.at(i)->Rebin(rebin);
        //hist.at(i)->GetXaxis()->SetRangeUser(xmin, xmax);
        //hist.at(i)->Smooth(1);

        // Iterative fitting
        int peakBin = hist.at(i)->GetMaximumBin();
        double peakValue = hist.at(i)->GetXaxis()->GetBinCenter(peakBin);
        double totalIntegral = hist.at(i)->Integral();
        double targetIntegral = 0.95 * totalIntegral;

        int binLeft = peakBin, binRight = peakBin;
        double cumulativeSum = hist.at(i)->GetBinContent(peakBin);

        while (cumulativeSum < targetIntegral) {
            if (binLeft > 1) cumulativeSum += hist.at(i)->GetBinContent(--binLeft);
            if (binRight < hist.at(i)->GetNbinsX()) cumulativeSum += hist.at(i)->GetBinContent(++binRight);
            if (binLeft == 1 && binRight == hist.at(i)->GetNbinsX()) break;
        }

        double fitMin = hist.at(i)->GetXaxis()->GetBinCenter(binLeft);
        double fitMax = hist.at(i)->GetXaxis()->GetBinCenter(binRight);
        fitMin=xmin; fitMax=xmax;

        //TF1* fitFunc = new TF1("CruijffFit", Cruijff, fitMin, fitMax, 6);
        Fit_func[i] = new TF1("CruijffFit", Cruijff, fitMin, fitMax, 6);
        Fit_func[i]->SetParameters(peakValue, 0.01, 0.01, 0.1, 0.1, hist.at(i)->GetMaximum());
        Fit_func[i]->SetParLimits(1, 0.00001, 0.1);
        Fit_func[i]->SetParLimits(2, 0.00001, 0.1);

        for (int bin = 1; bin <= hist.at(i)->GetNbinsX(); ++bin) {
            double content = hist.at(i)->GetBinContent(bin);
            double error = hist.at(i)->GetBinError(bin);

            if (content > 0 && error > 0) {
                double weight = content / error;  // Weight = 1 / relative error
                hist.at(i)->SetBinError(bin, 1.0 / weight);  // Inverse weight becomes the new error
            }
        }
        double old_mean=hist.at(i)->GetMean();
        for (int iter = 0; iter < 7; ++iter) {
            hist.at(i)->Fit(Fit_func[i], "R W");
            double mean = Fit_func[i]->GetParameter(0);
            double sigmaL = Fit_func[i]->GetParameter(1);
            double sigmaR = Fit_func[i]->GetParameter(2);
            if ((abs(old_mean-mean)/mean) < 0.0002) break;
            //fitMin = mean - 2 * sigmaL;
            //fitMax = mean + 2 * sigmaR;
            
            
        }

        /*for (int iter = 0; iter < 7; ++iter) {
            hist.at(i)->Fit(Fit_func[i], "R E");
            double mean = Fit_func[i]->GetParameter(0);
            double sigmaL = Fit_func[i]->GetParameter(1);
            double sigmaR = Fit_func[i]->GetParameter(2);
            //fitMin = mean - 2 * sigmaL;
            //fitMax = mean + 2 * sigmaR;
            //Fit_func[i]->SetRange(fitMin, fitMax);
        }*/

        hist.at(i)->GetXaxis()->SetRangeUser(fitMin, fitMax);
        //hist.at(i)->GetYaxis()->SetRangeUser(0, 1.4 * hist.at(0)->GetMaximum());
        hist.at(i)->GetYaxis()->SetRangeUser(0, 1.5 *MaxY);
        hist.at(i)->SetMarkerColor(line_color[i]);
        hist.at(i)->SetMarkerStyle(marker_style[i]);
        hist.at(i)->SetMarkerSize(2);
        Fit_func[i]->SetLineColor(line_color[i]);
        //if(i==0){ hist.at(i)->Draw("same");}
        //else hist.at(i)->Draw("same");
        //Fit_func[i]->Draw("SAME");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextColor(line_color[i]);
        latex.SetTextSize(0.035);
        latex.SetText(0.8, 0.9-(i*0.15), Form("#mu = %.4f", Fit_func[i]->GetParameter(0)));
        latexVec[i].push_back(latex);
        latex.SetText(0.8, 0.85-(i*0.15), Form("#sigma_{L} = %.4f", Fit_func[i]->GetParameter(1)));
        latexVec[i].push_back(latex);
        latex.SetText(0.8, 0.8-(i*0.15), Form("#sigma_{R} = %.4f", Fit_func[i]->GetParameter(2)));
        latexVec[i].push_back(latex);
        
        legend->AddEntry(hist.at(i), legend_texts[i].c_str(), "p");
    }
    for(int i = 0; i < (int)latexVec.size(); i++) {
    for(int j = 0; j < (int)latexVec[i].size(); j++) {
        latexVec[i][j].Draw();  // Draw each TLatex object
    }
        hist.at(i)->Draw("same");
        Fit_func[i]->Draw("same");
}
    legend->Draw();
    gPad->Modified();
    gPad->Update();
    cout<<Fit_func[0]->GetParameter(0)<<endl;
    cout<<Fit_func[1]->GetParameter(0)<<endl;

    if (save_canvas) {
        canvas_n1->SaveAs(Form("%s.png", tag_name));
        canvas_n1->SaveAs(Form("%s.pdf", tag_name));
    }
}

struct MixedData {
     string str1;
     string str4;
     string str2;
     int intData;
     double double1;
     double double2;
     double double3;
     double double4;
     string str3;
 };

void Multifile_fit(){

    char* hname = new char[200];
  
  char* hist_name = new char[200];
  
  char* title= new char[2000];
 
  char *leg_head = new char[200];
 
  int n=0;
  int n_files=1;
 
    //f[0] = new TFile("plot.root");
    
    vector<string> filetag=  {"DRN response","BDT response"};
    //vector<vector<string>> varName;
    //vector<vector<string>> legend_texts;
    //vector<string> xLabel;
    vector<string> loghist;
    vector<string> norm;
MixedData varName[]= {  // { names of plots,Title of plot, xlabel, rebin, ymin, ymax , xmin, xmax, legend}
{"DRN_corrbyGen","DRN_ErecoByEgen","E_{Reco}/E_{Gen}",200,0,1,0.95,1.05,""},
{"BDT_corrbyGen","BDT_ErecoByEgen","E_{Reco}/E_{Gen}",200,0,1,0.95,1.05,""},

};
 vector<string> GEN = {"Angle between gen photons"};

loghist = {"Leading A boost", "Subleading A boost"};
norm ={"ErecoByEgen","ErecobyEgen"};
  sprintf(hname,"temp.root");
  TFile* fout = new TFile(hname,"RECREATE");
 vector<string> File_list={"Plot_UND4.root","Plot_UND2.root","Plot_UND3.root","Plot_UND1.root"};
    n_files=File_list.size();
  for(int i=0; i<size(varName); i++)
    {      
      //for(int i_file=0; i_file<n_files;i_file++){
        
        int rebin = varName[i].intData; 
        string xLabel = varName[i].str2;
        double ymin = varName[i].double1;
        double ymax = varName[i].double2;
        double xmin = varName[i].double3;
        double xmax = varName[i].double4;
        //vector<string> legend_texts = varName[i].str3;
        vector<string> legend_texts={"Run2","Threshold varied","Noise varied","Both varied"};
        string Name = varName[i].str4;
        string  VarName = varName[i].str1;

        vector<TH1F*> hist_list;
        for (int j=0; j<File_list.size();j++){
          f[j] = new TFile(File_list[j].c_str());
 
	  sprintf(hist_name,"%s",VarName.c_str());
	  cout<<hist_name<<"\t"<<i<<"\t"<<j<<"\t"<<f[j]->GetName()<<endl;
          
	  TH1F* h_resp2 = (TH1F*)f[j]->Get(hist_name); // SR
	  h_resp2->GetXaxis()->SetTitle(xLabel.c_str());
	  cout<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin<<"\t"<<xmin<<"\t"<<xmax<< "\t" <<"File_name "<<File_list[j]<<endl;
	  
	  //h_resp2->Rebin(rebin);
	 
	  
	  //h_resp2= setMyRange(h_resp2,xmin,xmax);
	  //setLastBinAsOverFlow(h_resp2);
	  
	  
	  hist_list.push_back(h_resp2); 
        }
      string  Savename;
    int gen = count(GEN.begin(),GEN.end(),Name);
    int LOG = count(loghist.begin(), loghist.end(),Name);
    int NORM= count(norm.begin(), norm.end(),Name);
    if(gen){Savename = "GEN_";}
    else {Savename = "RECO_";}
         string Savename2=Savename+to_string(1000 + i)+"_" +Name + "_overlay";

generate_1Dplot(hist_list,Savename2.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,false,false,false,true,filetag[i].c_str(),legend_texts);
if(LOG && NORM){
         string Savename1=Savename +to_string(1000 + i) +"_"+Name + "_overlay_logy_norm";
generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,true,true,false,true,filetag[i].c_str(),legend_texts);
}

          else if(LOG && !NORM){
          string Savename1=Savename+to_string(1000 + i)+"_" +Name + "_overlay_logy";
generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,false,true,false,true,filetag[i].c_str(),legend_texts);
	  }
          else if(!LOG && NORM){
          string Savename1=Savename+to_string(1000 + i)+"_" +Name + "_overlay_norm"; 
generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,true,false,false,true,filetag[i].c_str(),legend_texts);
}
	         
 }
}