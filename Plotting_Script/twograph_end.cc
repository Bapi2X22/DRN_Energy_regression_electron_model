#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TMath.h>
#include <vector>

// Cruijff function definition
double Cruijff(double *x, double *par) {
    double mean = par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double alphaL = par[3];
    double alphaR = par[4];
    double norm = par[5];

    double dx = x[0] - mean;
    double sigma = (dx < 0) ? (2*sigmaL * sigmaL + alphaL * dx * dx) : (2*sigmaR * sigmaR + alphaR * dx * dx);
    return norm * TMath::Exp(- dx * dx / sigma);
    
};


    double xmin = 0.95;
    double xmax = 1.05;
void superimposePlotsWithFits() {
    // std::vector<std::string> fileNames = {"Plot_UND4.root", "Plot_UND2.root", "Plot_UND3.root", "Plot_UND1.root"};
    std::vector<std::string> fileNames = {"Plot_end_4.root"};
    std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kBlack};
    // std::vector<std::string> labels = {"None Varied", "Rechit threshold varied", "Noise threshold varied", "Both varied"};
    std::vector<std::string> labels = {"BDT", "DRN"};
    // std::vector<int> Entries = {802294, 803486, 803240, 803498};
    
    
    TMultiGraph *mg = new TMultiGraph();
    // TMultiGraph *mg_DRN = new TMultiGraph();
    // TLegend *legend_BDT = new TLegend(0.7, 0.7, 0.9, 0.9);
    // TLegend *legend_DRN = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.02);

    std::vector<TF1*> fitFuncs_BDT;
    std::vector<TF1*> fitFuncs_DRN;
    std::vector<std::vector<double>> fitParams_BDT;
    std::vector<std::vector<double>> fitParams_DRN;

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
        
        BDT_ratio->Rebin(200);
        DRN_ratio->Rebin(200);

        BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        
        // DRN_ratio->GetYaxis()->SetRangeUser(0, 1.5 *MaxY);
        // BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        // DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
        DRN_ratio->Scale(1.0 / DRN_ratio->Integral());

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


            graph_BDT->SetPoint(j - 1, x, y_BDT);
            graph_BDT->SetPointError(j - 1, 0, err_BDT);
    
            graph_DRN->SetPoint(j - 1, x, y_DRN);
            graph_DRN->SetPointError(j - 1, 0, err_DRN);

            BDT_ratio->SetBinError(j, err_BDT);
            DRN_ratio->SetBinError(j, err_DRN);
        }

        graph_BDT->SetMarkerStyle(20 + i);
        graph_BDT->SetMarkerColor(colors[i]);
        graph_BDT->SetLineColor(colors[i]);

        graph_DRN->SetMarkerStyle(21 + i);
        graph_DRN->SetMarkerColor(colors[i+1]);
        graph_DRN->SetLineColor(colors[i+1]);

        mg->Add(graph_BDT);
        mg->Add(graph_DRN);
        // legend_BDT->AddEntry(graph_BDT, labels[i].c_str(), "lp");
        // legend_DRN->AddEntry(graph_DRN, labels[i+1].c_str(), "lp");
        legend->AddEntry(graph_BDT, "BDT", "lp");
        legend->AddEntry(graph_DRN, "DRN", "lp");

                        // Initialize fit function
        TF1 *fitFunc_BDT = new TF1(Form("fit_BDT_%zu", i), Cruijff, xmin, xmax, 6);
        TF1 *fitFunc_DRN = new TF1(Form("fit_DRN_%zu", i), Cruijff, xmin , xmax, 6);
        fitFunc_BDT->SetParLimits(1, 0.00001, 0.1);
        fitFunc_BDT->SetParLimits(2, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(1, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(2, 0.00001, 0.1);
        
        fitFunc_BDT->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, BDT_ratio->GetMaximum());
        fitFunc_DRN->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, DRN_ratio->GetMaximum());

        double old_mean_BDT=BDT_ratio->GetMean();
        for (int iter = 0; iter < 40; ++iter) {
            BDT_ratio->Fit(fitFunc_BDT, "RW");
            double mean_BDT = fitFunc_BDT->GetParameter(0);
            double sigmaL_BDT = fitFunc_BDT->GetParameter(1);
            double sigmaR_BDT = fitFunc_BDT->GetParameter(2);
            if ((abs(old_mean_BDT-mean_BDT)/mean_BDT) < 0.0002) break;
        };
        double old_mean_DRN=DRN_ratio->GetMean();
        for (int iter_1 = 0; iter_1 < 1000; ++iter_1) {
            DRN_ratio->Fit(fitFunc_DRN, "RW");
            double mean_DRN = fitFunc_DRN->GetParameter(0);
            double sigmaL_DRN = fitFunc_DRN->GetParameter(1);
            double sigmaR_DRN = fitFunc_DRN->GetParameter(2);
            if ((abs(old_mean_DRN-mean_DRN)/mean_DRN) < 0.0002) break;
        };


        
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

        fitParams_BDT.push_back(params_BDT);
        fitParams_DRN.push_back(params_DRN);
        
        fitFunc_BDT->SetLineColor(colors[i]);
        fitFunc_BDT->SetLineStyle(i);
        fitFunc_DRN->SetLineColor(colors[i+1]);
        fitFunc_DRN->SetLineStyle(i);
        fitFuncs_BDT.push_back(fitFunc_BDT);
        fitFuncs_DRN.push_back(fitFunc_DRN);

        file->Close();
    }

    TCanvas *c_twoGraph = new TCanvas("c_twoGraph", "BDT vs DRN", 800, 600);
    c_twoGraph->cd();
    
    // Draw the BDT multigraph first
    mg->Draw("AP");
    mg->GetXaxis()->SetLimits(xmin, xmax);
    mg->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    mg->GetYaxis()->SetTitle("Normalized");
    mg->GetYaxis()->SetRangeUser(0, mg->GetYaxis()->GetXmax() * 0.14);
    mg->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->Draw("AP same");
    // mg_DRN->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->GetYaxis()->SetRangeUser(0, mg_DRN->GetYaxis()->GetXmax() * 0.14);

    // Draw the DRN multigraph on the same canvas
    
    // Add labels and legends
    TLatex *lateX = new TLatex();
    lateX->DrawTextNDC(0.19, 0.95, "BDT vs DRN response comparison (Endcap)");
    lateX->SetTextSize(0.03);
    lateX->DrawTextNDC(0.14, 0.55, "pfThresUL18_pedUL18");
    // legend_BDT->Draw();
    // legend_DRN->Draw("same");
    // TLegend *legend = new TLegend(0.65, 0.75, 0.85, 0.9);
    // legend->AddEntry(mg, "BDT", "P");
    // legend->AddEntry(mg_DRN, "DRN", "P");
    // legend->Draw();

    
    // Draw BDT fit functions
    int k = 0;
    for (auto fit : fitFuncs_BDT){
        fit->Draw("same");
        
        TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[k]);
        
        // Display fit parameters next to the fit line
        latex.DrawLatex(0.955, 0.13 + k * 0.03, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                                     fit->GetParameter(0),
                                                     fit->GetParameter(1),
                                                     fit->GetParameter(2)));
        k++;
    }
    
    // mg_DRN->Draw("P SAME");
    // mg_DRN->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->GetYaxis()->SetRangeUser(0, mg_DRN->GetYaxis()->GetXmax() * 0.14);
    // legend_DRN->Draw();
    // Draw DRN fit functions
    int l = 0; 
    for (auto fit : fitFuncs_DRN){ 
        fit->Draw("same");
        
        TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[l+1]);
        
        // Offset the label placement for DRN fits
        latex.DrawLatex(0.955, 0.13 + (l+1) * 0.01, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                                     fit->GetParameter(0),
                                                     fit->GetParameter(1),
                                                     fit->GetParameter(2)));
        l++;
    }
    legend->Draw();
    
    // Save the canvas
    c_twoGraph->SaveAs("BDT_vs_DRN_comparison(Endcap)_4.pdf");
    c_twoGraph->SaveAs("BDT_vs_DRN_comparison(Endcap)_4.png");

}
