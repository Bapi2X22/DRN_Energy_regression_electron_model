#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TMath.h>
#include <vector>
#include <fstream>  // For file writing

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


// Open output file
void superimposePlotsWithFits() {
    std::vector<std::string> fileNames = {"Photon_plot.root", "electron_plot.root"};
    // std::vector<std::string> fileNames = {"plot_pho.root", "plot_ele.root"};
    std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kBlack};
    std::vector<std::string> labels = {"Photon_model", "Electron_model"};
    // std::vector<int> Entries = {802294, 803486, 803240, 803498};
    
    
    TMultiGraph *mg_BDT = new TMultiGraph();
    TMultiGraph *mg_DRN = new TMultiGraph();
    TLegend *legend_BDT = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_DRN = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_BDT->SetTextSize(0.02);
    legend_DRN->SetTextSize(0.02);

    std::vector<TF1*> fitFuncs_BDT;
    std::vector<TF1*> fitFuncs_DRN;
    std::vector<std::vector<double>> fitParams_BDT;
    std::vector<std::vector<double>> fitParams_DRN;
    double chi2_ndf_BDT;
    double chi2_ndf_DRN;
    std::vector<double> chi2_ndf_vec;
    

    std::ofstream outFile("FittedParameters.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file!" << std::endl;
        return;
    }

    for (size_t i = 0; i < fileNames.size(); i++) {
        TFile *file = TFile::Open(fileNames[i].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open " << fileNames[i] << "!" << std::endl;
            continue;
        }

        TH1F *BDT_ratio = (TH1F*) file->Get("BDT_corrbyGen");
        TH1F *DRN_ratio = (TH1F*) file->Get("DRN_corrbyGen");

        if (!BDT_ratio || !DRN_ratio) {
            std::cerr << "Error: Missing histograms in " << fileNames[i] << "!" << std::endl;
            file->Close();
            continue;
        }

        BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
        DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
        
        BDT_ratio->Rebin(200);
        DRN_ratio->Rebin(200);


        
        // DRN_ratio->GetYaxis()->SetRangeUser(0, 1.5 *MaxY);
        // BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        // DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);


        // //////////// Iterative fitting///////////////////////////////
        // int peakBin_BDT = BDT_ratio->GetMaximumBin();
        // double peakValue_BDT = BDT_ratio->GetXaxis()->GetBinCenter(peakBin_BDT);


        // int peakBin_DRN = DRN_ratio->GetMaximumBin();
        // double peakValue_DRN = DRN_ratio->GetXaxis()->GetBinCenter(peakBin_DRN);


        int nBins = BDT_ratio->GetNbinsX();
        TGraphErrors *graph_BDT = new TGraphErrors(nBins);
        TGraphErrors *graph_DRN = new TGraphErrors(nBins);

        for (int j = 1; j <= nBins; j++) {
        double x = BDT_ratio->GetBinCenter(j);
        double y_BDT = BDT_ratio->GetBinContent(j);
        double y_DRN = DRN_ratio->GetBinContent(j);
        double bin_err_BDT = BDT_ratio->GetBinError(j);
        double bin_err_DRN = DRN_ratio->GetBinError(j);
        // double err_BDT = (y_BDT > 0) ? 1.0 / sqrt(Entries[i] * y_BDT) : 0.0;
        // double err_DRN = (y_DRN > 0) ? 1.0 / sqrt(Entries[i] * y_DRN) : 0.0;
        double err_BDT = (y_BDT > 0) ?  bin_err_BDT/y_BDT: 0.0;
        double err_DRN = (y_DRN > 0) ?  bin_err_DRN/y_DRN: 0.0;

        double width  = BDT_ratio->GetBinWidth(j);       // bin width
        double xErr   = width / 2.0;                    // symmetric error


        graph_BDT->SetPoint(j - 1, x, y_BDT);
        graph_BDT->SetPointError(j - 1, xErr, bin_err_BDT);

        graph_DRN->SetPoint(j - 1, x, y_DRN);
        graph_DRN->SetPointError(j - 1, xErr, bin_err_DRN);

        // BDT_ratio->SetBinError(j, err_BDT);
        // DRN_ratio->SetBinError(j, err_DRN);
    }

        graph_BDT->SetMarkerStyle(20 + i);
        graph_BDT->SetMarkerColor(colors[i]);
        graph_BDT->SetLineColor(colors[i]);
        graph_BDT->SetMarkerSize(0.5);
        

        graph_DRN->SetMarkerStyle(20 + i);
        graph_DRN->SetMarkerColor(colors[i]);
        graph_DRN->SetLineColor(colors[i]);
        graph_DRN->SetMarkerSize(0.5);

        mg_BDT->Add(graph_BDT);
        mg_DRN->Add(graph_DRN);
        legend_BDT->AddEntry(graph_BDT, labels[i].c_str(), "lp");
        legend_DRN->AddEntry(graph_DRN, labels[i].c_str(), "lp");

                        // Initialize fit function
        TF1 *fitFunc_BDT = new TF1(Form("fit_BDT_%zu", i), Cruijff, xmin, xmax, 6);
        TF1 *fitFunc_DRN = new TF1(Form("fit_DRN_%zu", i), Cruijff, xmin , xmax, 6);
        fitFunc_BDT->SetParLimits(1, 0.00001, 0.1);
        fitFunc_BDT->SetParLimits(2, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(1, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(2, 0.00001, 0.1);
        
        fitFunc_BDT->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, BDT_ratio->GetMaximum());
        fitFunc_DRN->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, DRN_ratio->GetMaximum());

        double old_mean_BDT = BDT_ratio->GetMean();
        double old_sigma_BDT = 0;
        
        for (int iter = 0; iter < 400; ++iter) {
            BDT_ratio->Fit(fitFunc_BDT, "RE");
        
            double mean_BDT = fitFunc_BDT->GetParameter(0);
            double sigmaL_BDT = fitFunc_BDT->GetParameter(1);
            double sigmaR_BDT = fitFunc_BDT->GetParameter(2);
            double sigma_BDT = (sigmaL_BDT+sigmaR_BDT)/2.0;
        
            bool mean_converged   = std::abs(old_mean_BDT - mean_BDT) / mean_BDT < 0.00002;
            bool sigma_converged = std::abs(old_sigma_BDT - sigma_BDT) / sigma_BDT < 0.00002;
        
            if (mean_converged && sigma_converged)
                break;
        
            old_mean_BDT = mean_BDT;
            old_sigma_BDT = sigma_BDT;
        }

        double old_mean_DRN = DRN_ratio->GetMean();
        double old_sigma_DRN = 0;
        
        for (int iter = 0; iter < 400; ++iter) {
            DRN_ratio->Fit(fitFunc_DRN, "RE");
        
            double mean_DRN = fitFunc_DRN->GetParameter(0);
            double sigmaL_DRN = fitFunc_DRN->GetParameter(1);
            double sigmaR_DRN = fitFunc_DRN->GetParameter(2);
            double sigma_DRN = (sigmaL_DRN+sigmaR_DRN)/2.0;
        
            bool mean_converged   = std::abs(old_mean_DRN - mean_DRN) / mean_DRN < 0.00002;
            bool sigma_converged = std::abs(old_sigma_DRN - sigma_DRN) / sigma_DRN < 0.00002;
        
            if (mean_converged && sigma_converged)
                break;
        
            old_mean_DRN = mean_DRN;
            old_sigma_DRN = sigma_DRN;
        }

        BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        
        // Iterative fitting process
        // double previousParams[6] = {1.0, 0.01, 0.01, 0.1, 0.1, BDT_ratio->GetMaximum()};
        // double accuracy = 1e-6;
        // int maxIterations = 100, iteration = 0;
        // double paramDifference = accuracy + 1.0;

        // while (paramDifference > accuracy && iteration < maxIterations) {
        //     iteration++;
        //     BDT_ratio->Fit(fitFunc_BDT, "RQ");
        //     DRN_ratio->Fit(fitFunc_DRN, "RQ");

        //     paramDifference = 0;
        //     for (int j = 0; j < 6; j++) {
        //         double diff_BDT = fabs(fitFunc_BDT->GetParameter(j) - previousParams[j]);
        //         double diff_DRN = fabs(fitFunc_DRN->GetParameter(j) - previousParams[j]);
        //         paramDifference += diff_BDT + diff_DRN;
        //         previousParams[j] = fitFunc_BDT->GetParameter(j);
        //     }
        // }

        // BDT_ratio->Fit(fitFunc_BDT, "RQ");
        // DRN_ratio->Fit(fitFunc_DRN, "RQ");
        std::vector<double> params_BDT;
        std::vector<double> params_DRN;

        for (int j = 0; j < 6; j++) {
            params_BDT.push_back(fitFunc_BDT->GetParameter(j));
            params_DRN.push_back(fitFunc_DRN->GetParameter(j));
        }

        double chi2_BDT = fitFunc_BDT->GetChisquare();
        int ndf_BDT = fitFunc_BDT->GetNDF();
        chi2_ndf_BDT = (ndf_BDT != 0) ? chi2_BDT / ndf_BDT : 0;  // Avoid division by zero

        double chi2_DRN = fitFunc_DRN->GetChisquare();
        int ndf_DRN = fitFunc_DRN->GetNDF();
        chi2_ndf_DRN = (ndf_DRN != 0) ? chi2_DRN / ndf_DRN : 0;  // Avoid division by zero
        chi2_ndf_vec.push_back(chi2_ndf_DRN);

        fitParams_BDT.push_back(params_BDT);
        fitParams_DRN.push_back(params_DRN);
        
        fitFunc_BDT->SetLineColor(colors[i]);
        fitFunc_DRN->SetLineColor(colors[i]);
        fitFuncs_BDT.push_back(fitFunc_BDT);
        fitFuncs_DRN.push_back(fitFunc_DRN);

        // outFile << "Dataset: " << labels[i] << "\n\n";

        //         // BDT Fit parameters
        // outFile << "BDT Fit Parameters:\n";
        // outFile << Form("  Mean: %.6f ± %.6f\n", fitFunc_BDT->GetParameter(0), fitFunc_BDT->GetParError(0));
        // outFile << Form("  SigmaL: %.6f ± %.6f\n", fitFunc_BDT->GetParameter(1), fitFunc_BDT->GetParError(1));
        // outFile << Form("  SigmaR: %.6f ± %.6f\n", fitFunc_BDT->GetParameter(2), fitFunc_BDT->GetParError(2));
        // outFile << Form("  AlphaL: %.6f ± %.6f\n", fitFunc_BDT->GetParameter(3), fitFunc_BDT->GetParError(3));
        // outFile << Form("  AlphaR: %.6f ± %.6f\n", fitFunc_BDT->GetParameter(4), fitFunc_BDT->GetParError(4));
        // outFile << Form("  Norm: %.6f ± %.6f\n", fitFunc_BDT->GetParameter(5), fitFunc_BDT->GetParError(5));

        // double chi2_BDT = fitFunc_BDT->GetChisquare();
        // int ndf_BDT = fitFunc_BDT->GetNDF();
        // double chi2_ndf_BDT = (ndf_BDT != 0) ? chi2_BDT / ndf_BDT : 0;  // Avoid division by zero

        // outFile << Form("  Chi²: %.6f\n", chi2_BDT);
        // outFile << Form("  NDF: %d\n", ndf_BDT);
        // outFile << Form("  Chi² / NDF: %.6f\n\n", chi2_ndf_BDT);

        // // DRN Fit parameters
        // outFile << "DRN Fit Parameters:\n";
        // outFile << Form("  Mean: %.6f ± %.6f\n", fitFunc_DRN->GetParameter(0), fitFunc_DRN->GetParError(0));
        // outFile << Form("  SigmaL: %.6f ± %.6f\n", fitFunc_DRN->GetParameter(1), fitFunc_DRN->GetParError(1));
        // outFile << Form("  SigmaR: %.6f ± %.6f\n", fitFunc_DRN->GetParameter(2), fitFunc_DRN->GetParError(2));
        // outFile << Form("  AlphaL: %.6f ± %.6f\n", fitFunc_DRN->GetParameter(3), fitFunc_DRN->GetParError(3));
        // outFile << Form("  AlphaR: %.6f ± %.6f\n", fitFunc_DRN->GetParameter(4), fitFunc_DRN->GetParError(4));
        // outFile << Form("  Norm: %.6f ± %.6f\n", fitFunc_DRN->GetParameter(5), fitFunc_DRN->GetParError(5));

        // double chi2_DRN = fitFunc_DRN->GetChisquare();
        // int ndf_DRN = fitFunc_DRN->GetNDF();
        // double chi2_ndf_DRN = (ndf_DRN != 0) ? chi2_DRN / ndf_DRN : 0;  // Avoid division by zero

        // outFile << Form("  Chi²: %.6f\n", chi2_DRN);
        // outFile << Form("  NDF: %d\n", ndf_DRN);
        // outFile << Form("  Chi² / NDF: %.6f\n\n", chi2_ndf_DRN);

        // outFile << "----------------------------------------\n";

 

        file->Close();
        
    }

    outFile.close();
    std::cout << "Fitted parameters, chi², NDF, and their ratio written to 'FittedParameters.txt'." << std::endl;



    // TCanvas *c_BDT = new TCanvas("c_BDT", "BDT Response Comparison", 800, 600);
    // mg_BDT->Draw("AP");

    // double yMax = 0;
    // TIter next(mg->GetListOfGraphs());
    // TGraph* g;
    // while ((g = (TGraph*)next())) {
    //     for (int i = 0; i < g->GetN(); ++i) {
    //         double x, y;
    //         g->GetPoint(i, x, y);
    //         if (y > yMax) yMax = y;
    //     }
    // }
    // mg->GetYaxis()->SetRangeUser(0, yMax * 1.25); // Add some margin (25%) on top
    // mg_BDT->GetXaxis()->SetLimits(xmin, xmax);
    // // mg_BDT->GetXaxis()->SetTitle("BDT Response");
    // mg_BDT->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    // mg_BDT->GetYaxis()->SetTitle("Normalized");
    // // mg_BDT->SetTitle("BDT response for 4 datasets");
    // TLatex *lateX = new TLatex();
    // lateX->DrawTextNDC(0.19, 0.95, "BDT response for 4 datasets");
    // // mg_BDT->GetYaxis()->SetRangeUser(0, mg_BDT->GetYaxis()->GetXmax() * 0.15);
    // legend_BDT->Draw();
    // int k = 0;
    // for (auto fit : fitFuncs_BDT){
    //     fit->Draw("same");
    //     TLatex latex;
    //     latex.SetTextSize(0.025);
    //     latex.SetTextColor(colors[k]);
    //     latex.DrawLatex(0.955, 0.12 + k * 0.01, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
    //                                                  fit->GetParameter(0),
    //                                                  fit->GetParameter(1),
    //                                                  fit->GetParameter(2)));
    //     k++;
    // };
    // c_BDT->SaveAs("BDT_Response_Comparison.pdf");
    // c_BDT->SaveAs("BDT_Response_Comparison.png");

    TCanvas *c_DRN = new TCanvas("c_DRN", "DRN Response Comparison", 800, 600);
    mg_DRN->Draw("AP");
    double yMax = 0;
    TIter next(mg_DRN->GetListOfGraphs());
    TGraph* g;
    while ((g = (TGraph*)next())) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (y > yMax) yMax = y;
        }
    }
    mg_DRN->GetYaxis()->SetRangeUser(0, yMax * 1.25); // Add some margin (25%) on top
    mg_DRN->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->GetXaxis()->SetTitle("E_corr/E_gen");
    // mg_DRN->SetTitle("DRN response for 4 datasets");
    TLatex *lateX = new TLatex();
    lateX->SetTextSize(0.034);
    lateX->DrawTextNDC(0.17, 0.95, "DRN_Response_Comparison_electron_vs_photon_model");
    lateX->SetTextSize(0.025);
    // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
    lateX->SetTextColor(colors[0]);
    lateX->DrawLatexNDC(0.17, 0.70, Form("#chi^{2}/NDF_photon = %.4f", chi2_ndf_vec[0]));
    lateX->SetTextColor(colors[1]);
    lateX->DrawLatexNDC(0.17, 0.65, Form("#chi^{2}/NDF_electron = %.4f", chi2_ndf_vec[1]));
    mg_DRN->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    mg_DRN->GetYaxis()->SetTitle("Normalized");
    // mg_DRN->GetYaxis()->SetRangeUser(0, mg_DRN->GetYaxis()->GetXmax() * 0.15);
    legend_DRN->Draw();
    int l = 0; 
    for (auto fit : fitFuncs_DRN){ 
        fit->Draw("same");
        TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[l]);
        double xTextPos = xmin + (xmax - xmin) * 0.05;  // 5% from left
        double yBase    = yMax * 1.2;                  // 95% of Y max (near top)
        double yStep    = yMax * 0.07;                  // spacing between lines
        
        // Display fit parameters next to the fit line
        latex.DrawLatex(xTextPos, yBase - l * yStep, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                              fit->GetParameter(0),
                                              fit->GetParameter(1),
                                              fit->GetParameter(2)));
        l++;
    };
    c_DRN->SaveAs("DRN_Response_Comparison_electron_vs_photon_model_HS.pdf");
    c_DRN->SaveAs("DRN_Response_Comparison_electron_vs_photon_model_HS.png");
}
