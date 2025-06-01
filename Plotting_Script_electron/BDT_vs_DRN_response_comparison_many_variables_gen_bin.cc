#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <iomanip>


using namespace std;

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

const double xmin = 0.95;
const double xmax = 1.05;


std::tuple<TGraphErrors*, TF1*, std::vector<double>, double> FitAndGraphHistogram(TH1F* hist, const std::string& dir, const std::string& name, int color, double xmin, double xmax) {
    // if (!hist) return std::make_tuple(nullptr, nullptr, std::vector<double>());

    hist->Scale(1.0 / hist->Integral());
    hist->Rebin(20);

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
    fitFunc->SetParLimits(1, 0.00001, 0.1);
    fitFunc->SetParLimits(2, 0.00001, 0.1);
    fitFunc->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist->GetMaximum());

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

    // ---------------------- DRAWING SECTION ----------------------
    TCanvas* c1 = new TCanvas(("c_" + name).c_str(), ("Fit: " + name).c_str(), 800, 600);
    c1->cd();

    double yMax = hist->GetMaximum();

    graph->GetYaxis()->SetRangeUser(0, yMax * 1.25);

    graph->Draw("AP");
    fitFunc->Draw("same");
    graph->GetXaxis()->SetLimits(xmin, xmax);
    graph->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    graph->GetYaxis()->SetTitle("Normalized");

    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextColor(color);

    double xText = xmin + (xmax - xmin) * 0.05;
    // double yMax = hist->GetMaximum();
    double yBase = yMax * 1.2;
    double yStep = yMax * 0.07;

    latex.DrawLatex(xText, yBase - 0 * yStep, Form("#mu = %.4f, #sigma_{L} = %.4f, #sigma_{R} = %.4f",
                                                    fitFunc->GetParameter(0),
                                                    fitFunc->GetParameter(1),
                                                    fitFunc->GetParameter(2)));

    latex.DrawLatex(xText, yBase - 1 * yStep, Form("#chi^{2}/NDF = %.4f", chi2_ndf));

    c1->SaveAs(Form("%sFit_%s.png", dir.c_str(), name.c_str()));
    c1->SaveAs(Form("%sFit_%s.pdf", dir.c_str(), name.c_str()));
        // ------------------------------------------------------------

    return std::make_tuple(graph, fitFunc, params, chi2_ndf);


}

template <typename T>

void ProcessGraphs(const std::string& fileName,
                   const std::string& tag, // "B" or "D"
                //    const std::vector<std::pair<int, int>>& bins,
                   const std::vector<std::pair<T, T>>& bins,
                   const std::string& directory,
                   int color,
                   const std::string& var,
                   std::vector<TGraphErrors*>& graphs_mu_out,
                   std::vector<TGraphErrors*>& graphs_sigma_out) {

    TFile* file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return;
    }

    int numBins = bins.size();
    TGraphErrors* g_mu = new TGraphErrors(numBins);
    TGraphErrors* g_sigma = new TGraphErrors(numBins);

    for (size_t i = 0; i < bins.size(); ++i) {
        const auto& bin = bins[i];
        auto binLow = bin.first;
        auto binHigh = bin.second;
        auto x = (binLow + binHigh) / 2.0;
        auto width = (binHigh - binLow) / 2.0;

        std::string histName = var + "_bin_" + tag + "_" + std::to_string(binLow) + "_" + std::to_string(binHigh);
        TH1F* hist = (TH1F*) file->Get(histName.c_str());
        if (!hist) continue;

        auto [graph, fit, params, chi2ndf] = FitAndGraphHistogram(hist, directory, histName, color, xmin, xmax);
        double mu = params[0];
        double sigma = (params[1] + params[2]) / 2.0;

        g_mu->SetPoint(i, x, mu);
        g_mu->SetPointError(i, width, 0);

        g_sigma->SetPoint(i, x, sigma / mu);
        g_sigma->SetPointError(i, width, 0);
    }

    graphs_mu_out.push_back(g_mu);
    graphs_sigma_out.push_back(g_sigma);
    file->Close();
}

// ProcessGraphs("your_file.root", "B", bins, kRed, graphs_mu_BDT, graphs_sigma_BDT);
// ProcessGraphs("your_file.root", "D", bins, kBlue, graphs_mu_DRN, graphs_sigma_DRN);
template <typename T>
void PlotMuAndSigmaFromFiles(
    const std::string& fileName,
    const std::vector<std::pair<T, T>>& bins,const std::string& object , const std::string& res_dir ,const std::string& dir , const std::string& label, const std::string& quantity, double muMin, double muMax, double muPMin, double muPMax, double sigMin, double sigMax, double sigPMin, double sigPMax, const std::string& xAxisLabel, const std::string& detector
) {
    std::vector<TGraphErrors*> graphs_mu_BDT, graphs_sigma_BDT;
    std::vector<TGraphErrors*> graphs_mu_DRN, graphs_sigma_DRN;

    // GenerateGraphsFromTwoFiles(
    //     fileName_p, fileName_e,
    //     bins,
    //     graphs_mu_BDT, graphs_sigma_BDT,
    //     graphs_mu_DRN_1, graphs_sigma_DRN_1,
    //     graphs_mu_DRN_2, graphs_sigma_DRN_2
    // );
    
    ProcessGraphs(fileName, "B", bins, dir ,kRed, quantity ,graphs_mu_BDT, graphs_sigma_BDT);
    ProcessGraphs(fileName, "D", bins, dir ,kBlue, quantity ,graphs_mu_DRN, graphs_sigma_DRN);
    // === Create canvas for mu with ratio ===
    TString canvasNameMu = "c_mu";
    TCanvas* c_mu = new TCanvas(canvasNameMu, canvasNameMu, 800, 800);

    // Create pads
    TPad* pad1_mu = new TPad("pad1_mu", "pad1_mu", 0, 0.3, 1, 1.0);
    TPad* pad2_mu = new TPad("pad2_mu", "pad2_mu", 0, 0.0, 1, 0.3);
    pad1_mu->SetBottomMargin(0.02);
    pad2_mu->SetTopMargin(0.02);
    pad2_mu->SetBottomMargin(0.3);
    pad1_mu->SetLeftMargin(0.15); // Default is ~0.1
    pad2_mu->SetLeftMargin(0.15); // Default is ~0.1
    pad1_mu->Draw();
    pad2_mu->Draw();

    // === Upper pad for mu ===
    pad1_mu->cd();
    pad1_mu->SetGrid();

    graphs_mu_BDT[0]->SetLineColor(kRed);
    graphs_mu_BDT[0]->SetMarkerStyle(20);
    graphs_mu_BDT[0]->SetTitle(Form("%s_ideal_IC_Sample", object.c_str()));
    graphs_mu_BDT[0]->SetMarkerColor(kRed);
    graphs_mu_BDT[0]->Draw("APL");
    graphs_mu_BDT[0]->GetYaxis()->SetRangeUser(muMin, muMax);
    graphs_mu_BDT[0]->GetYaxis()->SetTitle("Response(#mu)");
    graphs_mu_BDT[0]->GetXaxis()->SetLabelSize(0);

    graphs_mu_DRN[0]->SetMarkerColor(kBlue);
    graphs_mu_DRN[0]->SetLineColor(kBlue);
    graphs_mu_DRN[0]->SetMarkerStyle(21);
    graphs_mu_DRN[0]->Draw("PL SAME");

    TLatex* text_run3 = new TLatex();
    text_run3->SetNDC();                 // Use normalized coordinates
    text_run3->SetTextFont(62);          // Standard font
    text_run3->SetTextSize(0.035);       // Adjust size as needed
    text_run3->SetTextAlign(21);         // Horizontal center alignment
    text_run3->DrawLatex(0.5, 0.82, label.c_str());  // (x, y, text)


    // TLegend* leg_mu = new TLegend(0.6, 0.7, 0.88, 0.88);
    // leg_mu->AddEntry(graphs_mu_BDT[0], "BDT", "lp");
    // leg_mu->AddEntry(graphs_mu_DRN_1[0], "photon_DRN", "lp");
    // leg_mu->AddEntry(graphs_mu_DRN_2[0], "electron_DRN", "lp");
    // leg_mu->Draw();

    TLegend* leg_mu = new TLegend(0.65, 0.75, 0.85, 0.87); // Smaller box
    leg_mu->SetTextSize(0.025); // Smaller font size
    leg_mu->AddEntry(graphs_mu_BDT[0], "BDT", "lp");
    leg_mu->AddEntry(graphs_mu_DRN[0], "DRN", "lp");
    leg_mu->Draw();

    // === Lower pad for mu ratio ===
    pad2_mu->cd();
    int nPoints = graphs_mu_BDT[0]->GetN();
    TGraphErrors* ratio_BDT_mu = new TGraphErrors(nPoints);
    TGraphErrors* ratio_DRN_mu = new TGraphErrors(nPoints);

    for (int k = 0; k < nPoints; ++k) {
        double x, y_BDT, y_DRN;
        double err_BDT, err_DRN;

        graphs_mu_BDT[0]->GetPoint(k, x, y_BDT);
        graphs_mu_DRN[0]->GetPoint(k, x, y_DRN);

        err_BDT = graphs_mu_BDT[0]->GetErrorY(k);
        err_DRN = graphs_mu_DRN[0]->GetErrorY(k);

        ratio_DRN_mu->SetPoint(k, x, y_DRN / y_BDT);
        ratio_DRN_mu->SetPointError(k, 0, (y_DRN / y_BDT) * sqrt(pow(err_DRN / y_DRN, 2) + pow(err_BDT / y_BDT, 2)));
    }


    ratio_DRN_mu->SetMarkerColor(kBlue);
    ratio_DRN_mu->SetLineColor(kBlue);
    ratio_DRN_mu->SetMarkerStyle(21);
    ratio_DRN_mu->SetTitle(( ";" + xAxisLabel + ";Ratio to BDT" ).c_str());

    ratio_DRN_mu->Draw("APL");
    ratio_DRN_mu->GetYaxis()->SetTitleSize(0.07);
    ratio_DRN_mu->GetYaxis()->SetTitleOffset(0.8);
    ratio_DRN_mu->GetYaxis()->SetLabelSize(0.07);
    ratio_DRN_mu->GetXaxis()->SetLabelSize(0.08);
    ratio_DRN_mu->GetXaxis()->SetTitleSize(0.1);
    ratio_DRN_mu->GetXaxis()->SetTitleOffset(1.0);
    ratio_DRN_mu->SetMinimum(muPMin);
    ratio_DRN_mu->SetMaximum(muPMax);

    ratio_DRN_mu->Draw("PL SAME");

    gPad->Update();
    double xmin1 = gPad->GetUxmin();
    double xmax1 = gPad->GetUxmax();
    TLine* line1 = new TLine(xmin1, 1.0, xmax1, 1.0); // from x = xmin to x = xmax at y = 1
    line1->SetLineStyle(2); // 2 = dotted
    line1->SetLineWidth(1); // optional, line thickness
    line1->SetLineColor(kBlack); // or any other color
    line1->Draw("same");

    // Save the mu canvas
    c_mu->SaveAs(Form("%s%s_vs_mu_comparison_set_%s.png", res_dir.c_str(), quantity.c_str(), detector.c_str()));
    c_mu->SaveAs(Form("%s%s_vs_mu_comparison_set_%s.pdf", res_dir.c_str(), quantity.c_str(), detector.c_str()));


    // === Repeat everything similarly for sigma ===

    TString canvasNameSigma = "c_sigma";
    TCanvas* c_sigma = new TCanvas(canvasNameSigma, canvasNameSigma, 800, 800);
    TPad* pad1_sigma = new TPad("pad1_sigma", "pad1_sigma", 0, 0.3, 1, 1.0);
    TPad* pad2_sigma = new TPad("pad2_sigma", "pad2_sigma", 0, 0.0, 1, 0.3);
    pad1_sigma->SetBottomMargin(0.02);
    pad2_sigma->SetTopMargin(0.02);
    pad2_sigma->SetBottomMargin(0.3);
    pad1_sigma->SetLeftMargin(0.15); // Default is ~0.1
    pad2_sigma->SetLeftMargin(0.15); // Default is ~0.1
    pad1_sigma->Draw();
    pad2_sigma->Draw();

    pad1_sigma->cd();
    pad1_sigma->SetGrid();

    graphs_sigma_BDT[0]->SetMarkerColor(kRed);
    graphs_sigma_BDT[0]->SetLineColor(kRed);
    graphs_sigma_BDT[0]->SetMarkerStyle(20);
    graphs_sigma_BDT[0]->SetTitle(Form("%s_ideal_IC_Sample", object.c_str()));
    graphs_sigma_BDT[0]->Draw("APL");
    graphs_sigma_BDT[0]->GetYaxis()->SetRangeUser(sigMin, sigMax);
    graphs_sigma_BDT[0]->GetYaxis()->SetTitle("#sigma/#mu");
    graphs_sigma_BDT[0]->GetXaxis()->SetLabelSize(0);

    graphs_sigma_DRN[0]->SetMarkerColor(kBlue);
    graphs_sigma_DRN[0]->SetLineColor(kBlue);
    graphs_sigma_DRN[0]->SetMarkerStyle(21);
    graphs_sigma_DRN[0]->Draw("PL SAME");

    TLatex* text_run2 = new TLatex();
    text_run2->SetNDC();                 // Use normalized coordinates
    text_run2->SetTextFont(62);          // Standard font
    text_run2->SetTextSize(0.035);       // Adjust size as needed
    text_run2->SetTextAlign(21);         // Horizontal center alignment
    text_run2->DrawLatex(0.5, 0.82, label.c_str());  // (x, y, text)


    // TLegend* leg_sigma = new TLegend(0.6, 0.7, 0.88, 0.88);
    // leg_sigma->AddEntry(graphs_sigma_BDT[0], "BDT", "lp");
    // leg_sigma->AddEntry(graphs_sigma_DRN_1[0], "photon_DRN", "lp");
    // leg_sigma->AddEntry(graphs_sigma_DRN_2[0], "electron_DRN", "lp");
    // leg_sigma->Draw();

    // TLegend* leg_sigma = new TLegend(0.15, 0.75, 0.35, 0.87); // Left side, smaller box
    TLegend* leg_sigma = new TLegend(0.65, 0.75, 0.85, 0.87);
    leg_sigma->SetTextSize(0.025); // Smaller font size
    leg_sigma->AddEntry(graphs_sigma_BDT[0], "BDT", "lp");
    leg_sigma->AddEntry(graphs_sigma_DRN[0], "DRN", "lp");
    leg_sigma->Draw();


    pad2_sigma->cd();
    TGraphErrors* ratio_BDT_sigma = new TGraphErrors(nPoints);
    TGraphErrors* ratio_DRN_sigma = new TGraphErrors(nPoints);

    for (int m = 0; m < nPoints; ++m) {
        double X, Y_BDT, Y_DRN;
        double Err_BDT, Err_DRN;

        graphs_sigma_BDT[0]->GetPoint(m, X, Y_BDT);
        graphs_sigma_DRN[0]->GetPoint(m, X, Y_DRN);

        Err_BDT = graphs_sigma_BDT[0]->GetErrorY(m);
        Err_DRN = graphs_sigma_DRN[0]->GetErrorY(m);

        ratio_DRN_sigma->SetPoint(m, X, Y_DRN / Y_BDT);
        ratio_DRN_sigma->SetPointError(m, 0, (Y_DRN / Y_BDT) * sqrt(pow(Err_DRN / Y_DRN, 2) + pow(Err_BDT / Y_BDT, 2)));
    }

    ratio_DRN_sigma->SetMarkerColor(kBlue);
    ratio_DRN_sigma->SetLineColor(kBlue);
    ratio_DRN_sigma->SetMarkerStyle(21);
    ratio_DRN_sigma->SetTitle(( ";" + xAxisLabel + ";Ratio to BDT" ).c_str());

    ratio_DRN_sigma->Draw("APL");
    ratio_DRN_sigma->GetYaxis()->SetTitleSize(0.07);
    ratio_DRN_sigma->GetYaxis()->SetTitleOffset(0.6);
    ratio_DRN_sigma->GetYaxis()->SetLabelSize(0.07);
    ratio_DRN_sigma->GetXaxis()->SetLabelSize(0.08);
    ratio_DRN_sigma->GetXaxis()->SetTitleSize(0.1);
    ratio_DRN_sigma->GetXaxis()->SetTitleOffset(1.0);
    ratio_DRN_sigma->SetMinimum(sigPMin);
    ratio_DRN_sigma->SetMaximum(sigPMax);


    gPad->Update();
    double xmin2 = gPad->GetUxmin();
    double xmax2 = gPad->GetUxmax();
    TLine* line2 = new TLine(xmin2, 1.0, xmax2, 1.0); // from x = xmin to x = xmax at y = 1
    line2->SetLineStyle(2); // 2 = dotted
    line2->SetLineWidth(1); // optional, line thickness
    line2->SetLineColor(kBlack); // or any other color
    line2->Draw("same");

    // Save the sigma canvas
    c_sigma->SaveAs(Form("%s%s_vs_sigma_comparison_set_%s.png", res_dir.c_str(), quantity.c_str(), detector.c_str()));
    c_sigma->SaveAs(Form("%s%s_vs_sigma_comparison_set_%s.pdf", res_dir.c_str(), quantity.c_str(), detector.c_str()));
}

std::vector<std::pair<int, int>> pT_bins = {
    {20,40},{40, 60}, {60, 80}, {80, 100}, {100, 120}, {120, 140}, {140, 160},
        {160, 180}, {180, 200}, {200, 220}, {220, 240}, {240, 260},
        {260, 280}, {280, 300}
    };


std::vector<std::pair<double, double>> EB_eta_bins = {
    {0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4},
{0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8},
    {0.8, 0.9}, {0.9, 1.0}, {1.0, 1.1}, {1.1, 1.2}, {1.2, 1.3},{1.3, 1.442}};

std::vector<std::pair<double, double>> EE_eta_bins = {
    {1.566, 1.7}, {1.7, 1.8}, {1.8, 1.9}, {1.9, 2.0},
{2.0, 2.1}, {2.1, 2.2}, {2.2, 2.3}, {2.3, 2.4},
    {2.4, 2.5}};

const std::string& eta_dir_barrel_photon = "/eos/user/b/bbapi/www/Ideal_Photon/eta_response_barrel/";
const std::string& eta_dir_endcap_photon = "/eos/user/b/bbapi/www/Ideal_Photon/eta_response_endcap/";
const std::string& pt_dir_barrel_photon = "/eos/user/b/bbapi/www/Ideal_Photon/pT_response_barrel/";
const std::string& pt_dir_endcap_photon = "/eos/user/b/bbapi/www/Ideal_Photon/pT_response_endcap/";
const std::string& eta_dir_barrel_electron = "/eos/user/b/bbapi/www/Ideal_Electron/eta_response_barrel/";
const std::string& eta_dir_endcap_electron = "/eos/user/b/bbapi/www/Ideal_Electron/eta_response_endcap/";
const std::string& pt_dir_barrel_electron = "/eos/user/b/bbapi/www/Ideal_Electron/pT_response_barrel/";
const std::string& pt_dir_endcap_electron = "/eos/user/b/bbapi/www/Ideal_Electron/pT_response_endcap/";
const std::string& res_dir_photon = "/eos/user/b/bbapi/www/Ideal_Photon/Response_curve/";
const std::string& res_dir_electron = "/eos/user/b/bbapi/www/Ideal_Electron/Response_curve/";


PlotMuAndSigmaFromFiles<int>("Plot_ideal_photon_barrel.root", pT_bins,"Photon",res_dir_photon, pt_dir_barrel_photon,"Run2_ideal_IC (EB)", "pT", 0.99, 1.015, 0.998, 1.004, 0.0, 0.02, 0.9, 1.1, "pT[GeV/c]", "Barrel");
PlotMuAndSigmaFromFiles<int>("Plot_ideal_photon_endcap.root", pT_bins, "Photon",res_dir_photon, pt_dir_endcap_photon,"Run2_ideal_IC (EE)", "pT", 0.99, 1.015, 0.998, 1.004, 0.0, 0.02, 0.8, 1.1, "pT[GeV/c]", "Endcap");
PlotMuAndSigmaFromFiles<double>("Plot_ideal_photon_barrel.root", EB_eta_bins, "Photon",res_dir_photon , eta_dir_barrel_photon ,"Run2_ideal_IC (EB)", "eta", 0.99, 1.01, 0.990, 1.01, 0.0045, 0.0075, 0.8, 1.2, "|#eta|", "Barrel");
PlotMuAndSigmaFromFiles<double>("Plot_ideal_photon_endcap.root", EE_eta_bins, "Photon",res_dir_photon , eta_dir_endcap_photon ,"Run2_ideal_IC (EE)", "eta", 0.99, 1.01, 0.990, 1.01, 0.0045, 0.0075, 0.8, 1.2, "|#eta|", "Endcap");

PlotMuAndSigmaFromFiles<int>("Plot_electron_ideal_barrel.root", pT_bins, "Electron",res_dir_electron, pt_dir_barrel_electron,"Run2_ideal_IC (EB)", "pT", 0.99, 1.015, 0.998, 1.0075, 0.0, 0.02, 0.9, 1.1, "pT[GeV/c]", "Barrel");
PlotMuAndSigmaFromFiles<int>("Plot_electron_ideal_endcap.root", pT_bins, "Electron",res_dir_electron, pt_dir_endcap_electron,"Run2_ideal_IC (EE)", "pT", 0.99, 1.020, 0.998, 1.02, 0.0, 0.035, 0.8, 1.1, "pT[GeV/c]", "Endcap");
PlotMuAndSigmaFromFiles<double>("Plot_electron_ideal_barrel.root", EB_eta_bins, "Electron",res_dir_electron , eta_dir_barrel_electron ,"Run2_ideal_IC (EB)", "eta", 0.99, 1.01, 0.990, 1.01, 0.0045, 0.01, 0.8, 1.2, "|#eta|", "Barrel");
PlotMuAndSigmaFromFiles<double>("Plot_electron_ideal_endcap.root", EE_eta_bins, "Electron",res_dir_electron , eta_dir_endcap_electron ,"Run2_ideal_IC (EE)", "eta", 0.99, 1.01, 0.990, 1.01, 0.006, 0.015, 0.8, 1.2, "|#eta|", "Endcap");



