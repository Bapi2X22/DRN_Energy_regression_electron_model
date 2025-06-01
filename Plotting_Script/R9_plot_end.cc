#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>

void plot_R9() {
    // List of ROOT files to process
    // std::vector<std::string> fileNames = {"Plot_SE1.root","Plot_SE2.root","Plot_SE3.root","Plot_SE4.root"};  // Add more files if needed
    std::vector<std::string> fileNames = {"Plot_end_1.root", "Plot_end_2.root", "Plot_end_3.root", "Plot_end_4.root"};
    // std::vector<std::string> fileNames = {"Plot_E1.root"};
    std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2, kViolet, kYellow+2};

    std::vector<std::string> labels = {"Both varied", "Rechits threshold varied", "Noise Threshold varied", "None varied"};

    // Create canvas and legend
    TCanvas *c = new TCanvas("c", "num_rechits Distribution", 800, 600);
    c->SetLogy();
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.02);

    // To accumulate histograms for overlaying
    std::vector<TH1F*> histograms;

    for (size_t i = 0; i < fileNames.size(); i++) {
        TFile *file = TFile::Open(fileNames[i].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open " << fileNames[i] << "!" << std::endl;
            continue;
        }

        // Access the histogram directly
        TH1F *R9 = (TH1F*) file->Get("R9");
        if (!R9) {
            std::cerr << "Error: Missing histogram 'sum_En_rechits' in " << fileNames[i] << "!" << std::endl;
            file->Close();
            continue;
        }

        // Clone the histogram to avoid modifying the original
        TH1F *hist = (TH1F*) R9->Clone(Form("hist_%zu", i));
        hist->Rebin(10);
        hist->Scale(1.0 / hist->Integral());
        hist->SetDirectory(0);  // Prevents ROOT from auto-deleting the histogram when closing the file
        hist->SetStats(0);

        // Style the histogram
        // hist->SetLineColor(colors[i]);
        hist->SetLineColor(colors[i % colors.size()]);
        hist->SetLineWidth(2);

        // Add to the list and legend
        histograms.push_back(hist);
        legend->AddEntry(hist, labels[i].c_str(), "l");

        file->Close();
        delete file;
    }

    // Draw all histograms on the same canvas
    bool first = true;
    for (auto* hist : histograms) {
        if (first) {
            hist->Draw("HIST");  // Draw the first histogram
            first = false;
        } else {
            hist->Draw("HIST SAME");  // Overlay subsequent histograms
        }
    hist->GetXaxis()->SetTitle("Distribution of R9 (Endcap)");
    hist->GetYaxis()->SetTitle("Normalized");
    }


    // Draw legend
    legend->Draw();

    // Save the plot as an image
    c->SaveAs("R9_dist_Endcap.png");

    // Clean up
    for (auto* hist : histograms) {
        delete hist;
    }
    delete legend;
    delete c;
}