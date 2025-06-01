#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;


// Cruijff function definition
// double Cruijff(double *x, double *par) {
//     double mean = par[0];
//     double sigmaL = par[1];
//     double sigmaR = par[2];
//     double alphaL = par[3];
//     double alphaR = par[4];
//     double norm = par[5];

//     double dx = x[0] - mean;
//     double sigma = (dx < 0) ? (2 * sigmaL * sigmaL + alphaL * dx * dx)
//                              : (2 * sigmaR * sigmaR + alphaR * dx * dx);
//     return norm * TMath::Exp(-dx * dx / sigma);
// };

double Cruijff(double *x, double *par) {
    double mean   = par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double alphaL = par[3];
    double alphaR = par[4];
    double norm   = par[5];

    double dx = x[0] - mean;
    double sigma2;

    if (dx < 0)
        sigma2 = 2.0 * sigmaL * sigmaL + alphaL * dx * dx;
    else
        sigma2 = 2.0 * sigmaR * sigmaR + alphaR * dx * dx;

    if (sigma2 <= 0 || !std::isfinite(sigma2))
        return 0.0; // Protect against division by zero or NaNs

    double exponent = -dx * dx / sigma2;

    if (!std::isfinite(exponent)) 
        return 0.0; // Again, catch bad exponent values

    return norm * std::exp(exponent);
}

    double xmin = 0.95;
    double xmax = 1.05;
void plotSigmaOverMu() {

    std::vector<std::string> fileNames = {"Plot_pT_B_1.root", "Plot_pT_B_2.root", "Plot_pT_B_3.root", "Plot_pT_B_4.root"};
    std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kBlack};
    // std::vector<std::string> labels = {"BDT", "DRN"};
    std::vector<std::string> labels = {"Both varied", "Rechit threshold varied", "Noise threshold varied", "None varied"};
        std::vector<std::vector<double>> muvec_BDT;
        std::vector<std::vector<double>> sigvec_BDT;

        std::vector<std::vector<double>> muvec_DRN;
        std::vector<std::vector<double>> sigvec_DRN;
        std::vector<double> pT_values;


    for (int i = 0; i < fileNames.size(); i++) {
        TFile *file = TFile::Open(fileNames[i].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open " << fileNames[i] << "!" << std::endl;
            continue;
        }
    
    // TFile *file = TFile::Open("Plot_pT_B_1.root", "READ");
    // if (!file || file->IsZombie()) {
    //     std::cerr << "Error: Cannot open plot.root!" << std::endl;
    //     return;
    // }

// std::vector<std::pair<int, int>> pT_bins = {
//     {60, 80},
//     {80, 100}, {100, 120}, {120, 140}, {140, 160},
//     {160, 180}, {180, 200}, {200, 220}, {220, 240},
//     {240, 260}, {260, 280}, {280, 300}
// };

std::vector<std::pair<int, int>> pT_bins = {
    {40, 60}, {60, 80},
    {80, 100}, {100, 120}, {120, 140}, {140, 160},
    {160, 180}, {180, 200}, {200, 220}, {220, 240},
    {240, 260}, {260, 280}, {280, 300},
};

    
    std::vector<double> sigma_over_mu_BDT;
    std::vector<double> sigma_over_mu_DRN;
    std::vector<double> mu_hist_BDT;
    std::vector<double> mu_hist_DRN;

    for (const auto &bin : pT_bins) {
        std::string histName_BDT = "pT_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
        std::string histName_DRN = "pT_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);

        TH1F *hist_BDT = (TH1F*) file->Get(histName_BDT.c_str());
        TH1F *hist_DRN = (TH1F*) file->Get(histName_DRN.c_str());

        if (!hist_BDT || !hist_DRN) {
            std::cerr << "Missing histogram: " << histName_BDT << " or " << histName_DRN << std::endl;
            continue;
        }

        TMultiGraph *mg = new TMultiGraph();
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->SetTextSize(0.02);

        // std::vector<std::string> labels = {"BDT", "DRN"};

        // std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kBlack};

        std::vector<std::vector<double>> fitParams_BDT;
        std::vector<std::vector<double>> fitParams_DRN;

        hist_BDT->Rebin(200);
        hist_DRN->Rebin(200);

        hist_BDT->GetXaxis()->SetRangeUser(xmin, xmax);
        hist_DRN->GetXaxis()->SetRangeUser(xmin, xmax);
        
        // DRN_ratio->GetYaxis()->SetRangeUser(0, 1.5 *MaxY);
        // BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        // DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        hist_BDT->Scale(1.0 / hist_BDT->Integral());
        hist_DRN->Scale(1.0 / hist_DRN->Integral());

        int nBins = hist_BDT->GetNbinsX();
        TGraphErrors *graph_BDT = new TGraphErrors(nBins);
        TGraphErrors *graph_DRN = new TGraphErrors(nBins);

        for (int j = 1; j <= nBins; j++) {
        double x = hist_BDT->GetBinCenter(j);
        double y_BDT = hist_BDT->GetBinContent(j);
        double y_DRN = hist_DRN->GetBinContent(j);
        double bin_err_BDT = hist_BDT->GetBinError(j);
        double bin_err_DRN = hist_DRN->GetBinError(j);
        // double err_BDT = (y_BDT > 0) ? 1.0 / sqrt(Entries[i] * y_BDT) : 0.0;
        // double err_DRN = (y_DRN > 0) ? 1.0 / sqrt(Entries[i] * y_DRN) : 0.0;
        double err_BDT = (y_BDT > 0) ?  bin_err_BDT/y_BDT: 0.0;
        double err_DRN = (y_DRN > 0) ?  bin_err_DRN/y_DRN: 0.0;

        graph_BDT->SetPoint(j - 1, x, y_BDT);
        graph_BDT->SetPointError(j - 1, 0, err_BDT);

        graph_DRN->SetPoint(j - 1, x, y_DRN);
        graph_DRN->SetPointError(j - 1, 0, err_DRN);

        hist_BDT->SetBinError(j, err_BDT);
        hist_DRN->SetBinError(j, err_DRN);
    }

        graph_BDT->SetMarkerStyle(20);
        graph_BDT->SetMarkerColor(colors[0]);
        graph_BDT->SetLineColor(colors[0]);

        graph_DRN->SetMarkerStyle(21);
        graph_DRN->SetMarkerColor(colors[1]);
        graph_DRN->SetLineColor(colors[1]);

        mg->Add(graph_BDT);
        mg->Add(graph_DRN);
        // legend_BDT->AddEntry(graph_BDT, labels[i].c_str(), "lp");
        // legend_DRN->AddEntry(graph_DRN, labels[i+1].c_str(), "lp");
        legend->AddEntry(graph_BDT, "BDT", "lp");
        legend->AddEntry(graph_DRN, "DRN", "lp");

        TF1 *fitFunc_BDT = new TF1("fit_BDT", Cruijff, xmin, xmax, 6);
        TF1 *fitFunc_DRN = new TF1("fit_DRN", Cruijff, xmin , xmax, 6);
        fitFunc_BDT->SetParLimits(1, 0.00001, 0.1);
        fitFunc_BDT->SetParLimits(2, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(1, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(2, 0.00001, 0.1);
        
        fitFunc_BDT->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist_BDT->GetMaximum());
        fitFunc_DRN->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist_DRN->GetMaximum());

        double old_mean_BDT=hist_BDT->GetMean();
        for (int iter = 0; iter < 40; ++iter) {
            hist_BDT->Fit(fitFunc_BDT, "RW");
            double mean_BDT = fitFunc_BDT->GetParameter(0);
            double sigmaL_BDT = fitFunc_BDT->GetParameter(1);
            double sigmaR_BDT = fitFunc_BDT->GetParameter(2);
            double sigma_BDT = (sigmaL_BDT+sigmaR_BDT)/2.0;
            double sig_by_mu_BDT = sigma_BDT / mean_BDT;
            if ((abs(old_mean_BDT-mean_BDT)/mean_BDT) < 0.0002) break;
        };
        double old_mean_DRN=hist_DRN->GetMean();
        for (int iter_1 = 0; iter_1 < 40; ++iter_1) {
            hist_DRN->Fit(fitFunc_DRN, "RW");
            double mean_DRN = fitFunc_DRN->GetParameter(0);
            double sigmaL_DRN = fitFunc_DRN->GetParameter(1);
            double sigmaR_DRN = fitFunc_DRN->GetParameter(2);
            double sigma_DRN = (sigmaL_DRN+sigmaR_DRN)/2.0;
            double sig_by_mu_DRN = sigma_DRN / mean_DRN;
            if ((abs(old_mean_DRN-mean_DRN)/mean_DRN) < 0.0002) break;
        };
        if (i==0){
        pT_values.push_back((bin.first + bin.second) / 2.0);
        }
        // sigma_over_mu_BDT.push_back(sig_by_mu_BDT);
        // sigma_over_mu_DRN.push_back(sig_by_mu_DRN);
        double mean_DRN = fitFunc_DRN->GetParameter(0);
        double sigmaL_DRN = fitFunc_DRN->GetParameter(1);
        double sigmaR_DRN = fitFunc_DRN->GetParameter(2);
        double sigma_DRN = (sigmaL_DRN+sigmaR_DRN)/2.0;
        double sig_by_mu_DRN = sigma_DRN / mean_DRN;
        mu_hist_DRN.push_back(mean_DRN);
        sigma_over_mu_DRN.push_back(sig_by_mu_DRN);

        double mean_BDT = fitFunc_BDT->GetParameter(0);
        double sigmaL_BDT = fitFunc_BDT->GetParameter(1);
        double sigmaR_BDT = fitFunc_BDT->GetParameter(2);
        double sigma_BDT = (sigmaL_BDT+sigmaR_BDT)/2.0;
        double sig_by_mu_BDT = sigma_BDT / mean_BDT;
        mu_hist_BDT.push_back(mean_BDT);
        sigma_over_mu_BDT.push_back(sig_by_mu_BDT);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "pT bin [" << bin.first << ", " << bin.second << "]\n";
        std::cout << "  DRN: mu = " << mean_DRN
              << ", sigma_L = " << sigmaL_DRN
              << ", sigma_R = " << sigmaR_DRN
              << ", sigma/mu = " << sig_by_mu_DRN << '\n';
        std::cout << "  BDT: mu = " << mean_BDT
              << ", sigma_L = " << sigmaL_BDT
              << ", sigma_R = " << sigmaR_BDT
              << ", sigma/mu = " << sig_by_mu_BDT << "\n\n";



        std::vector<double> params_BDT;
        std::vector<double> params_DRN;



        muvec_DRN.push_back(mu_hist_DRN);
        sigvec_DRN.push_back(sigma_over_mu_DRN);

        muvec_BDT.push_back(mu_hist_BDT);
        sigvec_BDT.push_back(sigma_over_mu_BDT);
        
        for (int j = 0; j < 6; j++) {
            params_BDT.push_back(fitFunc_BDT->GetParameter(j));
            params_DRN.push_back(fitFunc_DRN->GetParameter(j));
        }

        fitParams_BDT.push_back(params_BDT);
        fitParams_DRN.push_back(params_DRN);
        
        fitFunc_BDT->SetLineColor(colors[0]);
        fitFunc_BDT->SetLineStyle(0);
        fitFunc_DRN->SetLineColor(colors[1]);
        fitFunc_DRN->SetLineStyle(0);
        // fitFuncs_BDT.push_back(fitFunc_BDT);
        // fitFuncs_DRN.push_back(fitFunc_DRN);

    TCanvas *c_twoGraph = new TCanvas("c_twoGraph", "BDT vs DRN", 800, 600);
    c_twoGraph->cd();
    
    // Draw the BDT multigraph first
    mg->Draw("AP");
    mg->GetXaxis()->SetLimits(xmin, xmax);
    mg->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    mg->GetYaxis()->SetTitle("Normalized");
    mg->GetYaxis()->SetRangeUser(0, mg->GetYaxis()->GetXmax() * 0.25);
    mg->GetXaxis()->SetLimits(xmin, xmax);
    
    // Add labels and legends
    TLatex *lateX = new TLatex();
    lateX->DrawTextNDC(0.19, 0.95, "BDT_vs_DRN(Barrel)_test");
    lateX->SetTextSize(0.03);
    lateX->DrawTextNDC(0.17, 0.55, "Test");
    
    // Draw BDT fit functions
        fitFunc_BDT->Draw("same");
        
        TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[0]);
        
        // Display fit parameters next to the fit line
        latex.DrawLatex(0.955, 0.235 + 0 * 0.02, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                                     fitFunc_BDT->GetParameter(0),
                                                     fitFunc_BDT->GetParameter(1),
                                                     fitFunc_BDT->GetParameter(2)));
        fitFunc_DRN->Draw("same");
        
        // TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[0+1]);
        
        // Offset the label placement for DRN fits
        latex.DrawLatex(0.955, 0.235 + (0+1) * 0.02, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                                     fitFunc_DRN->GetParameter(0),
                                                     fitFunc_DRN->GetParameter(1),
                                                     fitFunc_DRN->GetParameter(2)));
    legend->Draw();
    
    // Save the canvas
    // c_twoGraph->SaveAs("BDT_vs_DRN_Barrel_test.pdf");
    c_twoGraph->SaveAs(Form("/eos/home-b/bbapi/www/Energy_regression/response_test_%d_barr/BDT_vs_DRN_Barrel_test_%d.png",(i+1) ,(int)bin.first));


}


    file->Close();
}

    // Create graphs

    // TCanvas *c = new TCanvas("c", "Sigma/mu vs pT BDT", 800, 600);
    // int k = 0;
    // for (auto part : sigvec_BDT){
    //     TGraph *graph_sigma_BDT = new TGraph(pT_values.size(), pT_values.data(), part.data());
    //     if (k==0){
    //     graph_sigma_BDT->Draw("APL");
    //     } else {
    //     graph_sigma_BDT->Draw("same");
    //         };
    //     graph_sigma_BDT->SetMarkerStyle(20);
    //     graph_sigma_BDT->SetMarkerColor(colors[k]);
    //     if (k==0){
    // graph_sigma_BDT->SetTitle("#sigma/#mu vs pT(Barrel) BDT");
    // graph_sigma_BDT->GetXaxis()->SetTitle("pT [GeV]");
    // graph_sigma_BDT->GetYaxis()->SetTitle("#sigma/#mu");
    // graph_sigma_BDT->GetYaxis()->SetRangeUser(0, 0.015);

    // TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // // legend->AddEntry(graph_sigma_BDT, labels[k], "p");
    // // legend->AddEntry(graph_sigma_DRN, "DRN #sigma/#mu", "p");
    // // legend->Draw();
    //     }
    //     k++;
    // }

    // c->SaveAs("sigma_over_mu_vs_pT_barrel_BDT.png");
    // c->SaveAs("sigma_over_mu_vs_pT_barrel_BDT.pdf");

    // TCanvas *c1 = new TCanvas("c1", "Sigma/mu vs pT DRN", 800, 600);
    // // c->SetLogy();
    // int l = 0;
    // for (auto part : sigvec_DRN){
    //     TGraph *graph_sigma_DRN = new TGraph(pT_values.size(), pT_values.data(), part.data());
    //     if (l==0){
    //     graph_sigma_DRN->Draw("APL");
    //     } else {
    //     graph_sigma_DRN->Draw("same");
    //         };
    //     graph_sigma_DRN->SetMarkerStyle(20);
    //     graph_sigma_DRN->SetMarkerColor(colors[l]);
    //     if (l==0){
    // graph_sigma_DRN->SetTitle("#sigma/#mu vs pT(Barrel) DRN");
    // graph_sigma_DRN->GetXaxis()->SetTitle("pT [GeV]");
    // graph_sigma_DRN->GetYaxis()->SetTitle("#sigma/#mu");
    // graph_sigma_DRN->GetYaxis()->SetRangeUser(0, 0.015);
    // // TLatex *lateX = new TLatex();

    // TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // // legend1->AddEntry(graph_sigma_DRN, labels[l], "p");
    // // legend->AddEntry(graph_sigma_DRN, "DRN #sigma/#mu", "p");
    // // legend1->Draw();
    //     }
    //     l++;
    // }

    // c1->SaveAs("sigma_over_mu_vs_pT_barrel_DRN.png");
    // c1->SaveAs("sigma_over_mu_vs_pT_barrel_DRN.pdf");

    // TCanvas *c2 = new TCanvas("c2", "mu vs pT BDT", 800, 600);
    // gPad->SetLeftMargin(0.15); // Increase margin to create space


    // int m = 0;
    // for (auto part : muvec_BDT){
    //     TGraph *graph_mu_BDT = new TGraph(pT_values.size(), pT_values.data(), part.data());
    //     if (m==0){
    //     graph_mu_BDT->Draw("APL");
    //     } else {
    //     graph_mu_BDT->Draw("same");
    //         };
    //     graph_mu_BDT->SetMarkerStyle(20);
    //     graph_mu_BDT->SetMarkerColor(colors[m]);
    //     if (m==0){
    // graph_mu_BDT->SetTitle("##mu vs pT(Barrel) BDT");
    // graph_mu_BDT->GetXaxis()->SetTitle("pT [GeV]");
    // graph_mu_BDT->GetYaxis()->SetTitle("#sigma/#mu");
    // graph_mu_BDT->GetYaxis()->SetRangeUser(0.99, 1.01);
    // // TLatex *lateX = new TLatex();

    // TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // // legend2->AddEntry(graph_mu_BDT, labels[m], "p");
    // // legend->AddEntry(graph_sigma_DRN, "DRN #sigma/#mu", "p");
    // // legend2->Draw();
    //     }
    //     m++;
    // }
    // // c1->SetLogy();

    // c2->SaveAs("mu_vs_pT_barrel_BDT.pdf");
    // c2->SaveAs("mu_vs_pT_barrel_BDT.png");

    // TCanvas *c3 = new TCanvas("c3", "mu vs pT DRN", 800, 600);
    // gPad->SetLeftMargin(0.15); // Increase margin to create space

    // int n = 0;
    // for (auto part : muvec_DRN){
    //     TGraph *graph_mu_DRN = new TGraph(pT_values.size(), pT_values.data(), part.data());
    //     if (n==0){
    //     graph_mu_DRN->Draw("APL");
    //     } else {
    //     graph_mu_DRN->Draw("same");
    //         };
    //     graph_mu_DRN->SetMarkerStyle(20);
    //     graph_mu_DRN->SetMarkerColor(colors[n]);
    //     if (n==0){
    // graph_mu_DRN->SetTitle("##mu vs pT(Barrel) DRN");
    // graph_mu_DRN->GetXaxis()->SetTitle("pT [GeV]");
    // graph_mu_DRN->GetYaxis()->SetTitle("#sigma/#mu");
    // graph_mu_DRN->GetYaxis()->SetRangeUser(0.99, 1.01);
    // // TLatex *lateX = new TLatex();

    // TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // // legend3->AddEntry(graph_mu_DRN, labels[n], "p");
    // // legend->AddEntry(graph_sigma_DRN, "DRN #sigma/#mu", "p");
    // // legend3->Draw();
    //     }
    //     n++;
    // }

    // c3->SaveAs("mu_vs_pT_barrel_DRN.pdf");
    // c3->SaveAs("mu_vs_pT_barrel_DRN.png");

// ------------------ Canvas 0: Sigma/mu vs pT BDT ------------------
TCanvas *c = new TCanvas("c", "Sigma/mu vs pT BDT", 800, 600);
gPad->SetLeftMargin(0.15);
TMultiGraph *mg0 = new TMultiGraph();
std::vector<TGraph*> graphs_BDT;

int k = 0;
for (const auto& part : sigvec_BDT) {
    TGraph *g = new TGraph(pT_values.size(), pT_values.data(), part.data());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(colors[k]);
    // g->SetLineColor(colors[k]);
    mg0->Add(g, "AP");
    graphs_BDT.push_back(g);
    k++;
}

mg0->Draw("A");
mg0->SetTitle("#sigma/#mu vs pT (Barrel) BDT");
mg0->GetXaxis()->SetTitle("pT [GeV]");
mg0->GetYaxis()->SetTitle("#sigma/#mu");
mg0->GetYaxis()->SetRangeUser(0, 0.015);

TLegend *legend0 = new TLegend(0.7, 0.7, 0.9, 0.9);
for (int i = 0; i < 4; ++i) {
    legend0->AddEntry(mg0, labels[k].c_str(), "lp");
}
legend0->Draw();

c->SaveAs("sigma_over_mu_vs_pT_barrel_BDT.png");
c->SaveAs("sigma_over_mu_vs_pT_barrel_BDT.pdf");


// ------------------ Canvas 1: Sigma/mu vs pT DRN ------------------
TCanvas *c1 = new TCanvas("c1", "Sigma/mu vs pT DRN", 800, 600);
gPad->SetLeftMargin(0.15);
TMultiGraph *mg1 = new TMultiGraph();
std::vector<TGraph*> graphs_DRN;

int l = 0;
for (const auto& part : sigvec_DRN) {
    TGraph *g = new TGraph(pT_values.size(), pT_values.data(), part.data());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(colors[l]);
    // g->SetLineColor(colors[l]);
    mg1->Add(g, "AP");
    graphs_DRN.push_back(g);
    l++;
}

mg1->Draw("A");
mg1->SetTitle("#sigma/#mu vs pT (Barrel) DRN");
mg1->GetXaxis()->SetTitle("pT [GeV]");
mg1->GetYaxis()->SetTitle("#sigma/#mu");
mg1->GetYaxis()->SetRangeUser(0, 0.015);

TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
for (int i = 0; i < 4; ++i) {
    legend1->AddEntry(mg1, labels[k].c_str(), "lp");
}
legend1->Draw();

c1->SaveAs("sigma_over_mu_vs_pT_barrel_DRN.png");
c1->SaveAs("sigma_over_mu_vs_pT_barrel_DRN.pdf");


// ------------------ Canvas 2: mu vs pT BDT ------------------
TCanvas *c2 = new TCanvas("c2", "mu vs pT BDT", 800, 600);
gPad->SetLeftMargin(0.15);
TMultiGraph *mg2 = new TMultiGraph();
std::vector<TGraph*> muGraphs_BDT;

int m = 0;
for (const auto& part : muvec_BDT) {
    TGraph *g = new TGraph(pT_values.size(), pT_values.data(), part.data());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(colors[m]);
    // g->SetLineColor(colors[m]);
    mg2->Add(g, "AP");
    muGraphs_BDT.push_back(g);
    m++;
}

mg2->Draw("A");
mg2->SetTitle("#mu vs pT (Barrel) BDT");
mg2->GetXaxis()->SetTitle("pT [GeV]");
mg2->GetYaxis()->SetTitle("#mu");
mg2->GetYaxis()->SetRangeUser(0.99, 1.01);

TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
for (int i = 0; i < 4; ++i) {
    legend2->AddEntry(mg2, labels[k].c_str(), "lp");
}
legend2->Draw();

c2->SaveAs("mu_vs_pT_barrel_BDT.png");
c2->SaveAs("mu_vs_pT_barrel_BDT.pdf");


// ------------------ Canvas 3: mu vs pT DRN ------------------
TCanvas *c3 = new TCanvas("c3", "mu vs pT DRN", 800, 600);
gPad->SetLeftMargin(0.15);
TMultiGraph *mg3 = new TMultiGraph();
std::vector<TGraph*> muGraphs_DRN;

int n = 0;
for (const auto& part : muvec_DRN) {
    TGraph *g = new TGraph(pT_values.size(), pT_values.data(), part.data());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(colors[n]);
    // g->SetLineColor(colors[n]);
    mg3->Add(g, "AP");
    muGraphs_DRN.push_back(g);
    n++;
}

mg3->Draw("A");
mg3->SetTitle("#mu vs pT (Barrel) DRN");
mg3->GetXaxis()->SetTitle("pT [GeV]");
mg3->GetYaxis()->SetTitle("#mu");
mg3->GetYaxis()->SetRangeUser(0.99, 1.01);

TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
for (int i = 0; i < 4; ++i) {
    legend3->AddEntry(mg3, labels[k].c_str(), "lp");
}
legend3->Draw();

c3->SaveAs("mu_vs_pT_barrel_DRN.png");
c3->SaveAs("mu_vs_pT_barrel_DRN.pdf");

}
