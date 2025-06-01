#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

gROOT->SetBatch(kTRUE);

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

const double xmin = 80;
const double xmax = 100;


std::tuple<TGraphErrors*, TF1*, std::vector<double>, double> FitAndGraphHistogram(TH1F* hist, const std::string& name, int color, double xmin, double xmax) {
    // if (!hist) return std::make_tuple(nullptr, nullptr, std::vector<double>());

    // hist->Scale(1.0 / hist->Integral());
    hist->Rebin(250);

    int nBins = hist->GetNbinsX();
    TGraphErrors* graph = new TGraphErrors(nBins);

    for (int j = 1; j <= nBins; j++) {
        double x = hist->GetBinCenter(j);
        double y = hist->GetBinContent(j);
        double xErr = hist->GetBinWidth(j) / 2.0;
        double yErr = hist->GetBinError(j);
        graph->SetPoint(j - 1, x, y);
        graph->SetPointError(j - 1, xErr, yErr);
    }

    graph->SetMarkerStyle(20); 
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);


    TF1* fitFunc = new TF1(("fit_" + name).c_str(), Cruijff, xmin, xmax, 6);
    // fitFunc->SetParLimits(1, 0.00001, 0.1);
    // fitFunc->SetParLimits(2, 0.00001, 0.1);
    // fitFunc->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist->GetMaximum());

    fitFunc->SetParameters(91.0, 2.0, 2.0, 1.0, 1.0, hist->GetMaximum());
    fitFunc->SetParameters(90, 20, 20, 1, 1, hist->GetMaximum());

    fitFunc->SetParLimits(1, 0.000001, 5.0);   // sigmaL ~ [0.5, 5.0]
    fitFunc->SetParLimits(2, 0.000001, 5.0);   // sigmaR ~ [0.5, 5.0]

    double old_mean = hist->GetMean();
    double old_sigma = 0;
    for (int iter = 0; iter < 400; ++iter) {
        hist->Fit(fitFunc, "RE");
        double mean = fitFunc->GetParameter(0);
        double sigmaL = fitFunc->GetParameter(1);
        double sigmaR = fitFunc->GetParameter(2);
        double sigma = (sigmaL + sigmaR) / 2.0;

        bool mean_converged = std::abs(old_mean - mean) / mean < 0.00002;
        bool sigma_converged = std::abs(old_sigma - sigma) / sigma < 0.00002;

        if (mean_converged && sigma_converged) break;

        old_mean = mean;
        old_sigma = sigma;
    }

    hist->GetXaxis()->SetRangeUser(xmin, xmax);

    std::vector<double> params;
    for (int i = 0; i < 6; ++i)
        params.push_back(fitFunc->GetParameter(i));

    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    double chi2_ndf = (ndf != 0) ? chi2 / ndf : 0;

    fitFunc->SetLineColor(color);
    fitFunc->SetLineStyle(0);

    return std::make_tuple(graph, fitFunc, params, chi2_ndf);


}

    std::vector<std::string> fileName = {"plot_mee_40f_Photon.root", "plot_mee_40f.root"};

    TFile *file0 = TFile::Open(fileName[0].c_str(), "READ");
    TFile *file1 = TFile::Open(fileName[1].c_str(), "READ");

    TH1F *BDT_m_ee = (TH1F*) file0->Get("m_ee_BDT");
    TH1F *DRN_m_ee_1 = (TH1F*) file0->Get("m_ee_DRN");
    TH1F *DRN_m_ee_2 = (TH1F*) file1->Get("m_ee_DRN");

    std::string histName0 = "BDT_m_ee";
    std::string histName1 = "DRN_m_ee_1";
    std::string histName2 = "DRN_m_ee_2";

    auto [graph_BDT, fit_BDT, params_BDT, chi2ndf_BDT] = FitAndGraphHistogram(BDT_m_ee, histName0, kBlue, xmin, xmax);
    auto [graph_DRN_1, fit_DRN_1, params_DRN_1, chi2ndf_DRN_1] = FitAndGraphHistogram(DRN_m_ee_1, histName1, kMagenta, xmin, xmax);
    auto [graph_DRN_2, fit_DRN_2, params_DRN_2, chi2ndf_DRN_2] = FitAndGraphHistogram(DRN_m_ee_2, histName2, kRed, xmin, xmax);


    TCanvas* c1 = new TCanvas("c_threeGraph", "BDT vs electron_DRN vs photon_DRN", 800, 600);
    c1->cd();
    
    double yMax = std::max({BDT_m_ee->GetMaximum(), DRN_m_ee_1->GetMaximum(), DRN_m_ee_2->GetMaximum()});
    double yRange = yMax * 1.25;
    
    // Draw BDT
    graph_BDT->GetYaxis()->SetRangeUser(0, yRange);
    graph_BDT->GetXaxis()->SetLimits(xmin, xmax);
    graph_BDT->SetTitle("Invarient masss distribution (DYJetsToLL sample)");
    graph_BDT->GetXaxis()->SetTitle("m_ee [GeV/c^{2}]");
    graph_BDT->GetYaxis()->SetTitle("Events");
    graph_BDT->Draw("AP");
    fit_BDT->Draw("same");
    
    // Draw DRN 1
    graph_DRN_1->Draw("P SAME");
    fit_DRN_1->Draw("same");
    
    // Draw DRN 2
    graph_DRN_2->Draw("P SAME");
    fit_DRN_2->Draw("same");
    
    // Add legend
    // TLegend* legend = new TLegend(0.65, 0.65, 0.88, 0.88);
    TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.87);
    legend->SetTextSize(0.025); // Smaller font size
    legend->AddEntry(graph_BDT, "BDT", "p");
    legend->AddEntry(graph_DRN_1, "photon_DRN", "p");
    legend->AddEntry(graph_DRN_2, "electron_DRN", "p");
    legend->Draw();
    
    // Add text annotations (only for BDT here, but you can extend)
    TLatex latex;
    latex.SetTextSize(0.025);
    double xText = xmin + 0.05 * (xmax - xmin);
    double yStep = yMax * 0.07;
    double yBase = yMax * 1.2;

        // BDT text
    latex.SetTextColor(kBlue);
    latex.DrawLatex(xText, yBase - 0 * yStep, Form("#mu = %.4f, #sigma_{L} = %.4f, #sigma_{R} = %.4f",
        fit_BDT->GetParameter(0),
        fit_BDT->GetParameter(1),
        fit_BDT->GetParameter(2)));

    // DRN v1 text
    latex.SetTextColor(kMagenta);
    latex.DrawLatex(xText, yBase - 1 * yStep, Form("#mu = %.4f, #sigma_{L} = %.4f, #sigma_{R} = %.4f",
        fit_DRN_1->GetParameter(0),
        fit_DRN_1->GetParameter(1),
        fit_DRN_1->GetParameter(2)));

    // DRN v2 text
    latex.SetTextColor(kRed);  // slightly darker green for readability
    latex.DrawLatex(xText, yBase - 2 * yStep, Form("#mu = %.4f, #sigma_{L} = %.4f, #sigma_{R} = %.4f",
        fit_DRN_2->GetParameter(0),
        fit_DRN_2->GetParameter(1),
        fit_DRN_2->GetParameter(2)));

    latex.SetTextSize(0.025);
    // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
    latex.SetTextColor(kMagenta);
    latex.DrawLatexNDC(0.17, 0.60, Form("#chi^{2}/NDF_DRN_p = %.4f", chi2ndf_DRN_1));

    latex.SetTextSize(0.025);
    // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
    latex.SetTextColor(kRed);
    latex.DrawLatexNDC(0.17, 0.65, Form("#chi^{2}/NDF_DRN_e = %.4f", chi2ndf_DRN_2));

    latex.SetTextSize(0.025);
    // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
    latex.SetTextColor(kBlue);
    latex.DrawLatexNDC(0.17, 0.70, Form("#chi^{2}/NDF_BDT = %.4f", chi2ndf_BDT));


    c1->SaveAs("BDT_vs_eDRN_vspDRN_comparison_DY_mee_electron.pdf");
    c1->SaveAs("BDT_vs_eDRN_vspDRN_comparison_DY_mee_electron.png");