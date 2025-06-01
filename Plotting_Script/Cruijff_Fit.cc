#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>

// Cruijff function definition
double Cruijff(double *x, double *par) {
    double mean = par[0];   // Mean
    double sigmaL = par[1]; // Left sigma
    double sigmaR = par[2]; // Right sigma
    double alphaL = par[3]; // Left asymmetry
    double alphaR = par[4]; // Right asymmetry
    double norm = par[5];   // Normalization factor

    double dx = x[0] - mean;
    double sigma = (dx < 0) ? (2*sigmaL * sigmaL + alphaL * dx * dx) : (2*sigmaR * sigmaR + alphaR * dx * dx);
    return norm * TMath::Exp(-dx * dx / sigma);
}

void fitAndPlot(TH1F *hist, const char *title, const char *filename, Color_t color, double accuracy = 1e-4, int maxIterations = 100) {
    // Rebin the histogram for better statistics
    hist->Rebin(200);
    hist->Scale(1.0 / hist->Integral());

    
    // Restrict the range to 0.95 - 1.05
    hist->GetXaxis()->SetRangeUser(0.95, 1.05);
    
    // Define fit function
    TF1 *fitFunc = new TF1(Form("fit_%s", filename), Cruijff, 0.95, 1.05, 6);
    fitFunc->SetParNames("Mean", "SigmaL", "SigmaR", "AlphaL", "AlphaR", "Norm");
    fitFunc->SetParameters(1.0, 0.01, 0.01, 0.005, 0.005, hist->GetMaximum());
    
    // Create a new histogram for weighted fitting
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; i++) {
        double binContent = hist->GetBinContent(i);
        double binError = (binContent > 0) ? 1.0 / sqrt(802294*binContent) : 0.0;
        // Apply weight as the inverse square of the bin error
        hist->SetBinError(i, binError);
    }

    // Set the maximum number of iterations and accuracy for convergence
    double previousParams[6];
    for (int i = 0; i < 6; i++) previousParams[i] = fitFunc->GetParameter(i);

    int iteration = 0;
    double paramDifference = accuracy + 1.0; // Initially set to a value larger than accuracy

    while (paramDifference > accuracy && iteration < maxIterations) {
        iteration++;
        
        // Fit the histogram using weighted fitting
        hist->Fit(fitFunc, "R", "", 0.95, 1.05); // "R" for range, no "W" flag since we manually set errors

        // Calculate the difference between the current and previous parameters
        paramDifference = 0;
        for (int i = 0; i < 6; i++) {
            paramDifference += fabs(fitFunc->GetParameter(i) - previousParams[i]);
            previousParams[i] = fitFunc->GetParameter(i);
        }

        // Print iteration info for debugging
        if (iteration % 10 == 0) {
            cout << "Iteration " << iteration << ", Parameter difference: " << paramDifference << endl;
        }
    }

    // After fitting, display the number of iterations and final parameter values
    cout << "Fitting completed after " << iteration << " iterations with parameter difference: " << paramDifference << endl;

    // Convert histogram to points with error bars
    TGraphErrors *graph = new TGraphErrors(nBins);
    for (int i = 1; i <= nBins; i++) {
        double binCenter = hist->GetBinCenter(i);
        double binContent = hist->GetBinContent(i);
        double binError = hist->GetBinError(i);  // Use the updated error (weights)
        graph->SetPoint(i - 1, binCenter, binContent);
        graph->SetPointError(i - 1, 0, binError);
    }

    // Get min and max values for better axis scaling
    double minY = graph->GetYaxis()->GetXmin();
    double maxY = graph->GetYaxis()->GetXmax();
    graph->GetYaxis()->SetRangeUser(minY - 0.1 * fabs(minY), maxY + 0.1 * fabs(maxY));

    graph->GetXaxis()->SetLimits(0.95, 1.05);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(color);
    graph->GetYaxis()->SetRangeUser(0, 0.15 * maxY); // Adjust as needed

    // Create canvas
    TCanvas *c = new TCanvas(Form("c_%s", filename), title, 800, 600);
    graph->Draw("AP"); // 'A' to draw axis, 'P' for points
    fitFunc->Draw("same");

    // Extract fit parameters
    double mu = fitFunc->GetParameter(0);
    double sigmaL = fitFunc->GetParameter(1);
    double sigmaR = fitFunc->GetParameter(2);

    // Get the margins of the canvas (default margins)
    double leftMargin = c->GetLeftMargin();
    double rightMargin = c->GetRightMargin();
    double topMargin = c->GetTopMargin();
    
    // Draw the parameters (mu, sigmaL, sigmaR) on the canvas with four digits of precision
    TLatex *latex = new TLatex();
    latex->SetNDC(); // Use normalized coordinates
    latex->SetTextSize(0.035); // Font size
    latex->SetTextAlign(31);   // Align right (top right corner)
    latex->SetTextColor(kBlack); // Color of the text

    // Adjust the positions inside the canvas margins
    double xPos = 1.0 - rightMargin - 0.05; // Shift left by 5% from the right margin
    latex->DrawLatex(xPos, 1.0 - topMargin - 0.05, Form("#mu = %.4f", mu));        // Adjusted X and Y position
    latex->DrawLatex(xPos, 1.0 - topMargin - 0.10, Form("#sigma_{L} = %.4f", sigmaL)); // Adjusted Y position
    latex->DrawLatex(xPos, 1.0 - topMargin - 0.15, Form("#sigma_{R} = %.4f", sigmaR)); // Adjusted Y position
    
    // Save the plot as PNG and PDF
    c->SaveAs(Form("%s.png", filename));
    c->SaveAs(Form("%s.pdf", filename));
}




// Function to load histograms and call fit function
void runFits() {
    // Open ROOT file
    TFile *file = TFile::Open("Plot_UND1.root", "READ");
    if (!file || file->IsZombie()) {
        cout << "Error: Cannot open Plot_D1.root!" << endl;
        return;
    }

    // Retrieve histograms
    TH1F *DRN_ratio = (TH1F*) file->Get("DRN_corrbyGen");
    TH1F *BDT_ratio = (TH1F*) file->Get("BDT_corrbyGen");

    // Check if histograms exist
    if (!DRN_ratio || !BDT_ratio) {
        cout << "Error: One or both histograms not found!" << endl;
        file->Close();
        return;
    }

    // Fit and plot separately
    fitAndPlot(DRN_ratio, "DRN Fit", "DRN_Fit", kRed, 1e-5, 200);  // Example with 1e-5 accuracy and max 200 iterations
    fitAndPlot(BDT_ratio, "BDT Fit", "BDT_Fit", kBlue, 1e-5, 200);  // Example with 1e-5 accuracy and max 200 iterations


    // fitAndPlot(DRN_ratio, "DRN Fit", "DRN_Fit", kRed);
    // fitAndPlot(BDT_ratio, "BDT Fit", "BDT_Fit", kBlue);
    
    // Close file
    file->Close();
}


