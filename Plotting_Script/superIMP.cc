#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>

// Cruijff function definition
double Cruijff(double *x, double *par) {
    double mean = par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double alphaL = par[3];
    double alphaR = par[4];
    double norm = par[5];

    double dx = x[0] - mean;
    double sigma = (dx < 0) ? (sigmaL * sigmaL + alphaL * dx * dx) : (sigmaR * sigmaR + alphaR * dx * dx);
    return norm * TMath::Exp(-0.5 * dx * dx / sigma);
}

void fitAndPlotSuperimposed(const char* histName, const char* title, const char* filename, Color_t colors[4], const char* legendLabels[4], double accuracy = 1e-4, int maxIterations = 100) {
    const char* fileNames[4] = {"Plot_D1.root", "Plot_D2.root", "Plot_D3.root", "Plot_D4.root"};

    // Create a canvas
    TCanvas *c = new TCanvas(Form("c_%s", filename), title, 800, 600);
    TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.9); // Legend position

    bool firstPlot = true;

    for (int i = 0; i < 2; i++) {
        // Open the ROOT file
        TFile *file = TFile::Open(fileNames[i], "READ");
        if (!file || file->IsZombie()) {
            cout << "Error: Cannot open " << fileNames[i] << "!" << endl;
            continue;
        }

        // Retrieve histogram
        TH1F *hist = (TH1F*) file->Get(histName);
        if (!hist) {
            cout << "Error: Histogram " << histName << " not found in " << fileNames[i] << "!" << endl;
            file->Close();
            continue;
        }

        // Rebin and restrict range
        hist->Rebin(100);
        hist->GetXaxis()->SetRangeUser(0.95, 1.05);

        // Define fit function
        TF1 *fitFunc = new TF1(Form("fit_%s_%d", filename, i), Cruijff, 0.95, 1.05, 6);
        fitFunc->SetParNames("Mean", "SigmaL", "SigmaR", "AlphaL", "AlphaR", "Norm");
        fitFunc->SetParameters(1.0, 0.01, 0.01, 0.005, 0.005, hist->GetMaximum());

        // Weighted fitting
        int nBins = hist->GetNbinsX();
        for (int j = 1; j <= nBins; j++) {
            double binContent = hist->GetBinContent(j);
            double binError = (binContent > 0) ? 1.0 / sqrt(802294 * binContent) : 0.0;
            hist->SetBinError(j, binError);
        }

        // Iterative fitting
        double previousParams[6];
        for (int j = 0; j < 6; j++) previousParams[j] = fitFunc->GetParameter(j);

        int iteration = 0;
        double paramDifference = accuracy + 1.0;

        while (paramDifference > accuracy && iteration < maxIterations) {
            iteration++;
            hist->Fit(fitFunc, "R", "", 0.95, 1.05);

            paramDifference = 0;
            for (int j = 0; j < 6; j++) {
                paramDifference += fabs(fitFunc->GetParameter(j) - previousParams[j]);
                previousParams[j] = fitFunc->GetParameter(j);
            }
        }

        // Convert histogram to graph with error bars
        TGraphErrors *graph = new TGraphErrors(nBins);
        for (int j = 1; j <= nBins; j++) {
            double binCenter = hist->GetBinCenter(j);
            double binContent = hist->GetBinContent(j);
            double binError = hist->GetBinError(j);
            graph->SetPoint(j - 1, binCenter, binContent);
            graph->SetPointError(j - 1, 0, binError);
        }

        // double minY = graph->GetYaxis()->GetXmin();
        // double maxY = graph->GetYaxis()->GetXmax();
        // graph->GetYaxis()->SetLimits(-0.01, 0.14);
        // Set colors and markers
        double minY = graph->GetYaxis()->GetXmin();
        double maxY = graph->GetYaxis()->GetXmax();
        graph->GetYaxis()->SetRangeUser(minY - 0.1 * fabs(minY), maxY + 0.1 * fabs(maxY));
        graph->GetXaxis()->SetLimits(0.95, 1.05);
        graph->SetMarkerStyle(20 + i);
        graph->SetMarkerColor(colors[i]);
        graph->SetLineColor(colors[i]);
        graph->GetYaxis()->SetRangeUser(0, 0.1 * maxY); // Adjust as needed

        // Draw on canvas
        // c->cd();
        // if (firstPlot) {
        //     graph->Draw("AP");
        //     firstPlot = false;
        // } else {
        //     graph->Draw("P SAME");
        // }
        // fitFunc->Draw("same");

        // fitFunc->Draw("same");
        // delete fitFunc;
        graph->Draw("AP"); // 'A' to draw axis, 'P' for points
        fitFunc->Draw("same");
        // Add legend entry
        legend->AddEntry(graph, legendLabels[i], "lep");

        // Close file
        file->Close();
        // delete file;
    }

    // Draw legend and save plot
    c->cd();
    legend->Draw();
    c->SaveAs(Form("%s.png", filename));
    c->SaveAs(Form("%s.pdf", filename));
}

// Function to process all files and generate two superimposed plots
void runSuperimposedFits() {
    Color_t colors[4] = {kRed, kBlue, kGreen+2, kMagenta};
    const char* legendLabels[4] = {"Dataset 1", "Dataset 2", "Dataset 3", "Dataset 4"};

    // Superimpose BDT responses
    fitAndPlotSuperimposed("BDT_corrbyGen", "BDT Response Comparison", "BDT_Response_Comparison", colors, legendLabels, 1e-5, 200);

    // Superimpose DRN responses
    fitAndPlotSuperimposed("DRN_corrbyGen", "DRN Response Comparison", "DRN_Response_Comparison", colors, legendLabels, 1e-5, 200);
}