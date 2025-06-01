int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};                                                                              
int line_color[9] = {kBlue,kRed,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
vector<int> col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
void decorate(TH1F*,int);
void decorate(TH1F* hist,int i){
  hist->SetLineWidth(3);                                                                                                                        
}
void setLastBinAsOverFlow(TH1F*);
TH1F* setMyRange(TH1F*,double,double);

TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow                                                                                                                               
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;                                                                                                                     
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);                                                                                                                 
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
                                                                                                                         
  }
  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  cout<<xLow<<"\t"<<xHigh<<"\t"<<"set range"<<endl;
   return h1;
}

TH1F* DrawOverflow(TH1F*);
TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){

   UInt_t nx    = h->GetNbinsX()+1;
   Double_t *xbins= new Double_t[nx+1];
   for (UInt_t i=0;i<nx;i++)
     xbins[i]=h->GetBinLowEdge(i+1);
   xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
   char *tempName= new char[strlen(h->GetName())+10];
   sprintf(tempName,"%swtOverFlow",h->GetName());
   h->GetXaxis()->SetLimits(xmin,xrange);
   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
   htmp->GetXaxis()->SetRange(xmin,xrange);
   // Reset the axis labels
   htmp->SetXTitle(h->GetXaxis()->GetTitle());
   htmp->SetYTitle(h->GetYaxis()->GetTitle());
   // Fill the new hitogram including the extra bin for overflows
   for (UInt_t i=1; i<=nx; i++)
     htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
   // Fill the underflows
   htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   
   htmp->SetEntries(h->GetEntries());
   
   return htmp;
}
void setLastBinAsOverFlow(TH1F* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
  cout<<lastBinCt<<"\t"<<"Last bin values"<<endl;

}
void generate_1Dplot(vector<TH1F*> hist,char const *tag_name="",char const *xlabel="",char const *ylabel="",  int rebin=-1,double ymin=0,double ymax=0,int xmin=-1,int xmax=-1,char const *leg_head="",bool normalize=false, bool log_flag=false, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"})
{  
     TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,950,850);
       canvas_n1->Range(-60.25,-0.625,562.25,0.625);
       canvas_n1->SetFillColor(0);
       canvas_n1->SetBorderMode(0);
       canvas_n1->SetBorderSize(1);
       canvas_n1->SetLeftMargin(0.124);
       canvas_n1->SetRightMargin(0.035);
       canvas_n1->SetTopMargin(0.04);
       canvas_n1->SetBottomMargin(0.1);
       //THStack *hs_var=new THStack("var_Stack","");
       gStyle->SetOptStat(1111111);
       //   gStyle->SetOptStat(0);
       //double pvt_x_min = 0.6;
  double pvt_x_min = 0.75;
  double pvt_x_max = 0.99;
  double pvt_y_min = 0.9;
  //double pvt_dely = 0.18;
  double pvt_dely = 0.15;
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.21,0.82,0.65,0.95);  
  legend->SetTextSize(0.05);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(4);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title);
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.04;
  double xrange = xmax;
  
  vector<TH1F*> hist_list_temp;
  cout<<" hist.size() = "<<hist.size()<<endl;
  
  for(int i =0;i<(int)hist.size(); i ++) {
    
    if(normalize) {
      
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
      
    }
    
    hist.at(i)->SetLineWidth(line_width[i]);
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.11);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle("");
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.04);
    hist.at(i)->GetXaxis()->SetTitle(xlabel);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.040);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.2);
     decorate(hist.at(i),i);
    hist.at(i)->SetMarkerSize(0.8);
    hist.at(i)->SetMarkerStyle(20);
    hist.at(i)->SetMarkerColor(line_color[i]);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.3);
    hist.at(i)->GetXaxis()->SetLabelSize(0.04);

    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i].c_str(),"l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();

    

  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i = 0;i<(int)hist.size(); i++) {
    if(!normalize) {
    if (log_flag){ hist.at(i)->GetYaxis()->SetRangeUser(0.1,10*ymax);}
    else {
    hist.at(i)->GetYaxis()->SetRangeUser(0.0001,1.2*ymax);hist.at(i)->GetXaxis()->SetRangeUser(xmin,1.05*xmax);}
    }
    
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.0001,2.0);

      }
    cout<<"i"<<i<<endl;
    if(i==0) hist.at(i)->Draw("hist ");
    else hist.at(i)->Draw("hist sames");
	
  }
  legend->Draw();
  if(log_flag) {
      gPad->SetLogy();
    }
gPad->Update(); 
  TLatex* textOnTop = new TLatex();
  //new
    textOnTop->SetTextSize(0.054);
  //textOnTop->DrawLatexNDC(0.146,0.925,"CMS #it{#bf{Simulation Preliminary}}");
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.054);
  //textOnTop->DrawLatexNDC(0.72,0.925,en_lat);


  gPad->Modified();
                                                                                       
    
 char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);
     canvas_n1->SaveAs(canvas_name);   
     sprintf(canvas_name,"%s.pdf",tag_name);
    canvas_n1->SaveAs(canvas_name);  
  } 
}

struct MixedData {
    std::string str1;
    std::string str2;
    int intData;
    double double1;
    double double2;
    double double3;
    double double4;
    std::string str3;
};

const int nfiles=100;                                                                                                                                                             
TFile *f[nfiles];

void generate1Dplot()
{
  char* hname = new char[200];
  
  char* hist_name = new char[200];
  
  char* title= new char[2000];
 
  char *leg_head = new char[200];
 
  int n=0;
  int n_files=1;
 
    f[0] = new TFile("Plots_D1.root");
    
    vector<string> filetag=  {"Sample size:112K"};
MixedData varName[] = { // { Name of the plot , xLabel , rebin , ymin , ymax , xmin , xmax , Legend label }
{"DRN_corrbyGen","E_corr/E_gen",1,0,17,0.95,1.05,"BDT"} ,
{"BDT_corrbyGen","E_corr/E_gen",1,0,17,0.95,1.05,"DRN"} ,

};


  // vector<string>GEN ={"M_gen","A_gen_pT","A_gen_eta"}  ;                                                                                                               


//   bool flag=false;
 
//   sprintf(hname,"temp.root");
//   TFile* fout = new TFile(hname,"RECREATE");
 
//     n_files=1;  
    
//   for(int i_file=0; i_file<n_files;i_file++)
//     {      
     
//       for(int i_cut=0; i_cut<size(varName);i_cut++)
// 	{
//           int rebin = varName[i_cut].intData;
//           string xLabel = varName[i_cut].str2;
//           double ymin = varName[i_cut].double1;
//           double ymax = varName[i_cut].double2;
//           double xmin = varName[i_cut].double3;
//           double xmax = varName[i_cut].double4;
//           vector<string> legend_texts = {varName[i_cut].str3}; 
// 	  vector<TH1F*> hist_list;
	  
// 	  sprintf(hist_name,"%s",varName[i_cut].str1.c_str());

// 	  cout<<hist_name<<"\t"<<i_cut<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
          
// 	  TH1F* h_resp2 = (TH1F*)f[i_file]->Get(hist_name); // SR
// 	  h_resp2->GetXaxis()->SetTitle(xLabel.c_str());
// 	  cout<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin<<"\t"<<xmin<<"\t"<<xmax<<endl;
	  
// 	  h_resp2->Rebin(rebin);
	 
	  
// 	  h_resp2= setMyRange(h_resp2,xmin,xmax+0.01*xmax);
// 	  setLastBinAsOverFlow(h_resp2);
	  
	  
// 	  hist_list.push_back(h_resp2); 
	
	  	 
// 	  cout<<h_resp2->Integral()<<"\t"<<f[i_file]->GetName()<<endl;	 
// 	  cout<<" hist_list.size() "<<hist_list.size()<<"\t "<<endl;
//           std::ostringstream oss;
//           oss << std::setw(4) << std::setfill('0') << i_cut <<hist_name;
//           string hst_name = oss.str();

//           int gen = count(GEN.begin(),GEN.end(),varName[i_cut].str1);
// 	  int LOG = count(loghist.begin(), loghist.end(),varName[i_cut].str1);
//           if(gen){hst_name = "GEN_"+hst_name;}
//           else {hst_name = "RECO_" + hst_name;}
// 	  if(LOG){generate_1Dplot(hist_list,hst_name.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,false,true,false,true,filetag[i_file].c_str(),legend_texts);
// 	  }
// 	  else {
// generate_1Dplot(hist_list,hst_name.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,false,false,false,true,filetag[i_file].c_str(),legend_texts);
//           }
// 	}
//       //fout->Close();
      
//     }
    }
    
