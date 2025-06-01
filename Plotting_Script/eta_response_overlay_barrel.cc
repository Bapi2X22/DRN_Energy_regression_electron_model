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

void plotSigmaOverMuMulti() {
    std::vector<std::string> fileNames = {
        "Plot_ele_model_barrel_1.root",
        "Plot_ele_model_barrel_2.root",
        "Plot_ele_model_barrel_3.root",
        "Plot_ele_model_barrel_4.root"
    };

    std::vector<std::string> labels = {
        "Both varied", "Rechits varied", "Noise varied", "Run 2"
    };

    std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta+2};

    std::vector<std::pair<double, double>> EB_eta_bins = {
        {0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4},
    {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8},
        {0.8, 0.9}, {0.9, 1.0}, {1.0, 1.1}, {1.1, 1.2}, {1.2, 1.3},{1.3, 1.442}};

    std::vector<TGraphErrors*> graphs_sigma_BDT, graphs_sigma_DRN;
    std::vector<TGraphErrors*> graphs_mu_BDT, graphs_mu_DRN;

    for (int f = 0; f < fileNames.size(); ++f) {
        TFile* file = TFile::Open(fileNames[f].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open " << fileNames[f] << std::endl;
            continue;
        }

        std::vector<double> eta_values, sigma_over_mu_BDT, sigma_over_mu_DRN;
        std::vector<double> mu_hist_BDT, mu_hist_DRN;

        for (const auto& bin : EB_eta_bins) {
            std::string histName_BDT = "eta_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
            std::string histName_DRN = "eta_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);

            TH1F *hist_BDT = (TH1F*) file->Get(histName_BDT.c_str());
            TH1F *hist_DRN = (TH1F*) file->Get(histName_DRN.c_str());
            if (!hist_BDT || !hist_DRN) continue;

            TMultiGraph *mg = new TMultiGraph();
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->SetTextSize(0.02);

            std::vector<std::vector<double>> fitParams_BDT;
            std::vector<std::vector<double>> fitParams_DRN;

            hist_BDT->Scale(1.0 / hist_BDT->Integral());
            hist_DRN->Scale(1.0 / hist_DRN->Integral());

            hist_BDT->Rebin(50);
            hist_DRN->Rebin(50);
            
            // hist_BDT->GetXaxis()->SetRangeUser(xmin, xmax);
            // hist_DRN->GetXaxis()->SetRangeUser(xmin, xmax);

            // hist_BDT->Scale(1.0 / hist_BDT->Integral());
            // hist_DRN->Scale(1.0 / hist_DRN->Integral());

            int nBins = hist_BDT->GetNbinsX();
            TGraphErrors *graph_BDT = new TGraphErrors(nBins);
            TGraphErrors *graph_DRN = new TGraphErrors(nBins);
            
            for (int j = 1; j <= nBins; j++) {
                double x = hist_BDT->GetBinCenter(j);
                double y_BDT = hist_BDT->GetBinContent(j);
                double y_DRN = hist_DRN->GetBinContent(j);
                double width  = hist_BDT->GetBinWidth(j);       // bin width
                double xErr   = width / 2.0;                    // symmetric error
                double err_BDT = (y_BDT > 0) ? hist_BDT->GetBinError(j) / y_BDT : 0.0;
                double err_DRN = (y_DRN > 0) ? hist_DRN->GetBinError(j) / y_DRN : 0.0;
                double eBDT = hist_BDT->GetBinError(j);
                double eDRN = hist_DRN->GetBinError(j);
                graph_BDT->SetPoint(j - 1, x, y_BDT);
                graph_BDT->SetPointError(j - 1, xErr, eBDT);
        
                graph_DRN->SetPoint(j - 1, x, y_DRN);
                graph_DRN->SetPointError(j - 1, xErr, eDRN);
                // hist_BDT->SetBinError(j, err_BDT);
                // hist_DRN->SetBinError(j, err_DRN);
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
            TF1 *fitFunc_DRN = new TF1("fit_DRN", Cruijff, xmin, xmax, 6);
            fitFunc_BDT->SetParLimits(1, 0.00001, 0.1);
            fitFunc_BDT->SetParLimits(2, 0.00001, 0.1);
            fitFunc_DRN->SetParLimits(1, 0.00001, 0.1);
            fitFunc_DRN->SetParLimits(2, 0.00001, 0.1);

            fitFunc_BDT->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist_BDT->GetMaximum());
            fitFunc_DRN->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist_DRN->GetMaximum());

            // double old_mean_BDT = hist_BDT->GetMean();
            // for (int iter = 0; iter < 1000; ++iter) {
            //     hist_BDT->Fit(fitFunc_BDT, "RE");
            //     double mean = fitFunc_BDT->GetParameter(0);
            //     if (std::abs(old_mean_BDT - mean) / mean < 0.00002) break;
            // }

            // double old_mean_DRN = hist_DRN->GetMean();
            // for (int iter = 0; iter < 1000; ++iter) {
            //     hist_DRN->Fit(fitFunc_DRN, "RE");
            //     double mean = fitFunc_DRN->GetParameter(0);
            //     if (std::abs(old_mean_DRN - mean) / mean < 0.00002) break;
            // }

            double old_mean_BDT = hist_BDT->GetMean();
            double old_sigma_BDT = 0;
            
            for (int iter = 0; iter < 400; ++iter) {
                hist_BDT->Fit(fitFunc_BDT, "RE");
            
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

            double old_mean_DRN = hist_DRN->GetMean();
            double old_sigma_DRN = 0;
            
            for (int iter = 0; iter < 400; ++iter) {
                hist_DRN->Fit(fitFunc_DRN, "RE");
            
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


            hist_BDT->GetXaxis()->SetRangeUser(xmin, xmax);
            hist_DRN->GetXaxis()->SetRangeUser(xmin, xmax);
            double eta_center = (bin.first + bin.second) / 2.0;

            double mu_BDT = fitFunc_BDT->GetParameter(0);
            double sig_BDT = (fitFunc_BDT->GetParameter(1) + fitFunc_BDT->GetParameter(2)) / 2.0;
            mu_hist_BDT.push_back(mu_BDT);
            sigma_over_mu_BDT.push_back(sig_BDT / mu_BDT);

            double mu_DRN = fitFunc_DRN->GetParameter(0);
            double sig_DRN = (fitFunc_DRN->GetParameter(1) + fitFunc_DRN->GetParameter(2)) / 2.0;
            mu_hist_DRN.push_back(mu_DRN);
            sigma_over_mu_DRN.push_back(sig_DRN / mu_DRN);

            eta_values.push_back(eta_center);

            // std::cout << std::fixed << std::setprecision(4);
            // std::cout << "pT bin [" << bin.first << ", " << bin.second << "]\n";
            // std::cout << "  DRN: mu = " << mean_DRN
            //       << ", sigma_L = " << sigmaL_DRN
            //       << ", sigma_R = " << sigmaR_DRN
            //       << ", sigma/mu = " << sig_by_mu_DRN << '\n';
            // std::cout << "  BDT: mu = " << mean_BDT
            //       << ", sigma_L = " << sigmaL_BDT
            //       << ", sigma_R = " << sigmaR_BDT
            //       << ", sigma/mu = " << sig_by_mu_BDT << "\n\n";
    
    
    
            std::vector<double> params_BDT;
            std::vector<double> params_DRN;

            for (int j = 0; j < 6; j++) {
                params_BDT.push_back(fitFunc_BDT->GetParameter(j));
                params_DRN.push_back(fitFunc_DRN->GetParameter(j));
            }

            double chi2_BDT = fitFunc_BDT->GetChisquare();
            int ndf_BDT = fitFunc_BDT->GetNDF();
            double chi2_ndf_BDT = (ndf_BDT != 0) ? chi2_BDT / ndf_BDT : 0;  // Avoid division by zero

            double chi2_DRN = fitFunc_DRN->GetChisquare();
            int ndf_DRN = fitFunc_DRN->GetNDF();
            double chi2_ndf_DRN = (ndf_DRN != 0) ? chi2_DRN / ndf_DRN : 0;  // Avoid division by zero
    
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
        double yMax = 0;
        TIter next(mg->GetListOfGraphs());
        TGraph* g;
        while ((g = (TGraph*)next())) {
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                if (y > yMax) yMax = y;
            }
        }
        mg->GetYaxis()->SetRangeUser(0, yMax * 1.25); // Add some margin (25%) on top

        mg->Draw("AP");
        mg->GetXaxis()->SetLimits(xmin, xmax);
        mg->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
        mg->GetYaxis()->SetTitle("Normalized");
        // mg->GetYaxis()->SetRangeUser(0, mg->GetYaxis()->GetXmax() * 0.25);
        mg->GetXaxis()->SetLimits(xmin, xmax);
        
        // Add labels and legends
        TLatex *lateX = new TLatex();
        lateX->DrawTextNDC(0.19, 0.95, "BDT_vs_DRN(Barrel)_test");
        // Add labels and legends
        // TLatex *lateX = new TLatex();
        lateX->DrawTextNDC(0.19, 0.95, "BDT_vs_DRN(Barrel)_test");
        lateX->SetTextSize(0.025);
        // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
        lateX->SetTextColor(colors[1]);
        lateX->DrawLatexNDC(0.17, 0.65, Form("#chi^{2}/NDF_DRN = %.4f", chi2_ndf_DRN));

        lateX->SetTextSize(0.025);
        // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
        lateX->SetTextColor(colors[0]);
        lateX->DrawLatexNDC(0.17, 0.70, Form("#chi^{2}/NDF_BDT = %.4f", chi2_ndf_BDT));

        
        // Draw BDT fit functions
            fitFunc_BDT->Draw("same");
            
            TLatex latex;
            latex.SetTextSize(0.025);
            latex.SetTextColor(colors[0]);

            double xTextPos = xmin + (xmax - xmin) * 0.05;  // 5% from left
            double yBase    = yMax * 1.2;                  // 95% of Y max (near top)
            double yStep    = yMax * 0.07;                  // spacing between lines
            
            // Display fit parameters next to the fit line
            latex.DrawLatex(xTextPos, yBase - 0 * yStep, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                                  fitFunc_BDT->GetParameter(0),
                                                  fitFunc_BDT->GetParameter(1),
                                                  fitFunc_BDT->GetParameter(2)));
            fitFunc_DRN->Draw("same");
            
            // TLatex latex;
            latex.SetTextSize(0.025);
            latex.SetTextColor(colors[0+1]);
            
            // Offset the label placement for DRN fits
            latex.DrawLatex(xTextPos, yBase - 1 * yStep, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                                  fitFunc_DRN->GetParameter(0),
                                                  fitFunc_DRN->GetParameter(1),
                                                  fitFunc_DRN->GetParameter(2)));
        legend->Draw();
        
        // Save the canvas
        // c_twoGraph->SaveAs("BDT_vs_DRN_Barrel_test.pdf");
        // c_twoGraph->SaveAs(Form("/eos/home-b/bbapi/www/Energy_regression/response_test_%d_bar_new/BDT_vs_DRN_Barrel_test_"+std::to_string(bin.first) +"_"+std::to_string(bin.second)+".png",(f+1));
        //     }

        c_twoGraph->SaveAs(Form("/eos/home-b/bbapi/www/Energy_regression/eta_response_barrel/response_test_%d_barr/BDT_vs_DRN_Barrel_test_%.4f_%.4f.png", (f+1), (double)bin.first, (double)bin.second));

        }

        file->Close();
        
        int nPoints = eta_values.size(); // number of pT bins or eta bins
        TGraphErrors* Graph_mu_BDT = new TGraphErrors(nPoints);
        TGraphErrors* Graph_mu_DRN = new TGraphErrors(nPoints);
        TGraphErrors* Graph_sigma_BDT = new TGraphErrors(nPoints);
        TGraphErrors* Graph_sigma_DRN = new TGraphErrors(nPoints);
        
        for (int j = 0; j < nPoints; ++j) {
            // double x = (EB_eta_bins[j] + bin_edges[j+1]) / 2.0;     // bin center
            // double width = (bin_edges[j+1] - bin_edges[j]) / 2.0; // bin half-width
            double bin_low = EB_eta_bins[j].first;
            double bin_high = EB_eta_bins[j].second;
        
            double x = (bin_low + bin_high) / 2.0;         // bin center
            double width = (bin_high - bin_low) / 2.0;     // half bin width

            double y_BDT = mu_hist_BDT[j];
            double y_DRN = mu_hist_DRN[j];
        
            Graph_mu_BDT->SetPoint(j, x, y_BDT);
            Graph_mu_BDT->SetPointError(j, width, 0.0);  // x-error = bin width
        
            Graph_mu_DRN->SetPoint(j, x, y_DRN);
            Graph_mu_DRN->SetPointError(j, width, 0.0);
        }

        for (int j = 0; j < nPoints; ++j) {
            // double x = (EB_eta_bins[j] + bin_edges[j+1]) / 2.0;     // bin center
            // double width = (bin_edges[j+1] - bin_edges[j]) / 2.0; // bin half-width
            double bin_low = EB_eta_bins[j].first;
            double bin_high = EB_eta_bins[j].second;
        
            double x = (bin_low + bin_high) / 2.0;         // bin center
            double width = (bin_high - bin_low) / 2.0;     // half bin width

            double y_BDT = sigma_over_mu_BDT[j];
            double y_DRN = sigma_over_mu_DRN[j];
        
            Graph_sigma_BDT->SetPoint(j, x, y_BDT);
            Graph_sigma_BDT->SetPointError(j, width, 0.0);  // x-error = bin width
        
            Graph_sigma_DRN->SetPoint(j, x, y_DRN);
            Graph_sigma_DRN->SetPointError(j, width, 0.0);
        }


        // graphs_sigma_BDT.push_back(new TGraph(eta_values.size(), eta_values.data(), sigma_over_mu_BDT.data()));
        // graphs_sigma_DRN.push_back(new TGraph(eta_values.size(), eta_values.data(), sigma_over_mu_DRN.data()));
        // graphs_mu_BDT.push_back(new TGraph(eta_values.size(), eta_values.data(), mu_hist_BDT.data()));
        // graphs_mu_DRN.push_back(new TGraph(eta_values.size(), eta_values.data(), mu_hist_DRN.data()));
        graphs_sigma_BDT.push_back(Graph_sigma_BDT);
        graphs_sigma_DRN.push_back(Graph_sigma_DRN);
        graphs_mu_BDT.push_back(Graph_mu_BDT);
        graphs_mu_DRN.push_back(Graph_mu_DRN);
        
    }

    auto plotGraphs = [&](std::vector<TGraphErrors*> &graphs, const char* title, const char* yAxis, const char* filename, bool isMu) {
    TCanvas *c = new TCanvas(filename, title, 800, 800);

    // Upper pad for the main plot
    TPad *pad1 = new TPad("pad1", "Main", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    // TLegend *leg = new TLegend(0.12, 0.7, 0.35, 0.88);  // Top left
    TLegend *leg;
    if (std::string(yAxis).find("#sigma/#mu") != std::string::npos)
        leg = new TLegend(0.12, 0.7, 0.35, 0.88);  // Top left
    else
        leg = new TLegend(0.65, 0.7, 0.88, 0.88);  // Top right

    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetMarkerStyle(20 + i);
        graphs[i]->SetTitle(title);
        // graphs[i]->GetXaxis()->SetTitle("pT [GeV]");
        graphs[i]->GetXaxis()->SetLabelSize(0); // Hide tick labels
        graphs[i]->GetXaxis()->SetTitleSize(0);
        graphs[i]->GetYaxis()->SetTitle(yAxis);

        // Y-axis range
        if (std::string(yAxis).find("#sigma/#mu") != std::string::npos)
            graphs[i]->GetYaxis()->SetRangeUser(0, 0.014);
        else
            graphs[i]->GetYaxis()->SetRangeUser(0.99, 1.01);

        if (i == 0)
            graphs[i]->Draw("AP");
        else
            graphs[i]->Draw("P same");

        leg->AddEntry(graphs[i], labels[i].c_str(), "p");
    }
    leg->Draw();

    // Switch to canvas for ratio plot
    c->cd();
    TPad *pad2 = new TPad("pad2", "Ratio", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();

    std::vector<TGraph*> ratio_graphs;

    for (size_t i = 0; i < (graphs.size()-1); ++i) {
        std::vector<double> x_vals, ratio_vals;

        for (int j = 0; j < graphs[3]->GetN(); ++j) {
            double x, y0, y;
            graphs[3]->GetPoint(j, x, y0);
            graphs[i]->GetPoint(j, x, y);
            if (y0 != 0)
                ratio_vals.push_back(y / y0);
            else
                ratio_vals.push_back(0.0);
            x_vals.push_back(x);
        }

        TGraph *ratio = new TGraph(x_vals.size(), x_vals.data(), ratio_vals.data());
        ratio->SetLineColor(colors[i]);
        ratio->SetMarkerColor(colors[i]);
        ratio->SetMarkerStyle(20 + i);
        ratio_graphs.push_back(ratio);
    }

    double xmin = graphs[0]->GetXaxis()->GetXmin();
    double xmax = graphs[0]->GetXaxis()->GetXmax();

    // Dummy frame for axis
    TH1F *frame = pad2->DrawFrame(0, 0.9, 1.2, 1.5);   // Adjust as needed
    // TH1F *frame = nullptr;
    if (isMu) {
        TH1F *frame = pad2->DrawFrame(xmin, 0.995, xmax, 1.005);
    } else {
        TH1F *frame = pad2->DrawFrame(xmin, 0.9, xmax, 1.1);
    }
    frame->GetXaxis()->SetTitle("|#eta|");
    frame->GetYaxis()->SetTitle("Ratio w.r.t. Run2 sample");
    frame->GetXaxis()->SetTitleSize(0.10);
    frame->GetXaxis()->SetLabelSize(0.10);
    frame->GetYaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitleOffset(0.7);

    for (auto &g : ratio_graphs)
        g->Draw("PL same");

    c->SaveAs((std::string(filename) + ".pdf").c_str());
    c->SaveAs((std::string(filename) + ".png").c_str());
};

plotGraphs(graphs_sigma_BDT, "#sigma/#mu vs |#eta|_barrel (BDT)", "#sigma/#mu", "sigma_over_mu_vs_eta_BDT_barrel_ele_M", false);
plotGraphs(graphs_sigma_DRN, "#sigma/#mu vs |#eta|_barrel (DRN)", "#sigma/#mu", "sigma_over_mu_vs_eta_DRN_barrel_ele_M", false);
plotGraphs(graphs_mu_BDT, "#mu vs |#eta|_barrel (BDT)", "Response", "mu_vs_eta_BDT_barrel_ele_M", true);
plotGraphs(graphs_mu_DRN, "#mu vs |#eta|_barrel (DRN)", "Response", "mu_vs_eta_DRN_barrel_ele_M", true);

}
