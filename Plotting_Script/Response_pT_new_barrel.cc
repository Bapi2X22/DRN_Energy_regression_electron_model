#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <iostream>

// Cruijff function definition
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
    TFile *file = TFile::Open("Plot_pT_B_1.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open plot.root!" << std::endl;
        return;
    }

std::vector<std::pair<int, int>> pT_bins = {
    {40, 60}, {60, 80},
    {80, 100}, {100, 120}, {120, 140}, {140, 160},
    {160, 180}, {180, 200}, {200, 220}, {220, 240},
    {240, 260}, {260, 280}, {280, 300}
};

    std::vector<double> pT_values;
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

        hist_BDT->SetBinError(j, err_BDT);
        hist_DRN->SetBinError(j, err_DRN);
    }

        TF1 *fitFunc_BDT = new TF1("fit_BDT", Cruijff, xmin, xmax, 6);
        TF1 *fitFunc_DRN = new TF1("fit_DRN", Cruijff, xmin , xmax, 6);
        fitFunc_BDT->SetParLimits(1, 0.00001, 0.1);
        fitFunc_BDT->SetParLimits(2, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(1, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(2, 0.00001, 0.1);
        
        fitFunc_BDT->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist_BDT->GetMaximum());
        fitFunc_DRN->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist_DRN->GetMaximum());

        double old_mean_BDT=hist_BDT->GetMean();
        for (int iter = 0; iter < 1000; ++iter) {
            hist_BDT->Fit(fitFunc_BDT, "RW");
            double mean_BDT = fitFunc_BDT->GetParameter(0);
            double sigmaL_BDT = fitFunc_BDT->GetParameter(1);
            double sigmaR_BDT = fitFunc_BDT->GetParameter(2);
            double sigma_BDT = (sigmaL_BDT+sigmaR_BDT)/2.0;
            double sig_by_mu_BDT = sigma_BDT / mean_BDT;
            if ((abs(old_mean_BDT-mean_BDT)/mean_BDT) < 0.00002) break;
        };
        double old_mean_DRN=hist_DRN->GetMean();
        for (int iter_1 = 0; iter_1 < 1000; ++iter_1) {
            hist_DRN->Fit(fitFunc_DRN, "RW");
            double mean_DRN = fitFunc_DRN->GetParameter(0);
            double sigmaL_DRN = fitFunc_DRN->GetParameter(1);
            double sigmaR_DRN = fitFunc_DRN->GetParameter(2);
            double sigma_DRN = (sigmaL_DRN+sigmaR_DRN)/2.0;
            double sig_by_mu_DRN = sigma_DRN / mean_DRN;
            if ((abs(old_mean_DRN-mean_DRN)/mean_DRN) < 0.00002) break;
        };
        
        pT_values.push_back((bin.first + bin.second) / 2.0);
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
    }

    file->Close();

    // Create graphs
    TGraph *graph_sigma_BDT = new TGraph(pT_values.size(), pT_values.data(), sigma_over_mu_BDT.data());
    TGraph *graph_sigma_DRN = new TGraph(pT_values.size(), pT_values.data(), sigma_over_mu_DRN.data());

    TGraph *graph_mu_BDT = new TGraph(pT_values.size(), pT_values.data(), mu_hist_BDT.data());
    TGraph *graph_mu_DRN = new TGraph(pT_values.size(), pT_values.data(), mu_hist_DRN.data());

    graph_sigma_BDT->SetMarkerStyle(20);
    graph_sigma_BDT->SetMarkerColor(kRed);
    graph_sigma_DRN->SetMarkerStyle(21);
    graph_sigma_DRN->SetMarkerColor(kBlue);

    graph_mu_BDT->SetMarkerStyle(20);
    graph_mu_BDT->SetMarkerColor(kRed);
    graph_mu_DRN->SetMarkerStyle(21);
    graph_mu_DRN->SetMarkerColor(kBlue);

    TCanvas *c = new TCanvas("c", "Sigma/mu vs pT", 800, 600);
    // c->SetLogy();
    graph_sigma_BDT->Draw("APL");
    graph_sigma_DRN->Draw("PL same");

    graph_sigma_BDT->SetTitle("#sigma/#mu vs pT(Barrel)");
    graph_sigma_BDT->GetXaxis()->SetTitle("pT [GeV]");
    graph_sigma_BDT->GetYaxis()->SetTitle("#sigma/#mu");
    graph_sigma_BDT->GetYaxis()->SetRangeUser(0, 0.05);

    TLatex *lateX = new TLatex();
    lateX->SetTextSize(0.03);
    lateX->DrawTextNDC(0.45, 0.85, "Both_varied");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(graph_sigma_BDT, "BDT #sigma/#mu", "p");
    legend->AddEntry(graph_sigma_DRN, "DRN #sigma/#mu", "p");
    legend->Draw();

    c->SaveAs("sigma_over_mu_vs_pT_1_barrel.png");
    c->SaveAs("sigma_over_mu_vs_pT_1_barrel.pdf");

    TCanvas *c1 = new TCanvas("c1", "mu vs pT", 800, 600);
    gPad->SetLeftMargin(0.15); // Increase margin to create space

    // c1->SetLogy();
    graph_mu_BDT->Draw("APL");
    graph_mu_DRN->Draw("PL same");

    graph_mu_BDT->SetTitle("#mu vs pT(Barrel)");
    graph_mu_BDT->GetXaxis()->SetTitle("pT [GeV]");
    graph_mu_BDT->GetYaxis()->SetTitle("E_{corr}/E_{gen}");
    graph_mu_BDT->GetYaxis()->SetRangeUser(0.99, 1.01);

    // TLatex *lateX = new TLatex();
    lateX->SetTextSize(0.03);
    lateX->DrawTextNDC(0.45, 0.85, "Both_varied");

    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend1->AddEntry(graph_mu_BDT, "BDT #mu", "p");
    legend1->AddEntry(graph_mu_DRN, "DRN #mu", "p");
    legend1->Draw();

    c1->SaveAs("mu_vs_pT_1_barrel.pdf");
    c1->SaveAs("mu_vs_pT_1_barrel.png");

}
