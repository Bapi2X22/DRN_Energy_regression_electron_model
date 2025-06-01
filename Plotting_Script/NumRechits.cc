#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>

void plot_num_rechits() {
    // List of ROOT files to process
    std::vector<std::string> fileNames = {"Plot_all_1.root","Plot_all_2.root","Plot_all_3.root","Plot_all_4.root"};  // Add more files if needed
    std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2, kViolet, kYellow+2};

    std::vector<std::string> labels = {"Both varied", "Rechits threshold varied", "Noise Threshold varied", "None varied"};

    // Create canvas and legend
    TCanvas *c = new TCanvas("c", "num_rechits Distribution", 800, 600);
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
        TH1F *Num_Rechits = (TH1F*) file->Get("num_rechits");
        if (!Num_Rechits) {
            std::cerr << "Error: Missing histogram 'num_rechits' in " << fileNames[i] << "!" << std::endl;
            file->Close();
            continue;
        }

        // Clone the histogram to avoid modifying the original
        TH1F *hist = (TH1F*) Num_Rechits->Clone(Form("hist_%zu", i));
        hist->Rebin(10);
        hist->GetXaxis()->SetRangeUser(0, 120);
        hist->Scale(1.0 / hist->Integral());
        hist->SetDirectory(0);  // Prevents ROOT from auto-deleting the histogram when closing the file
        hist->SetStats(0);
        hist->SetTitle("");

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
    hist->GetXaxis()->SetTitle("Number of Rechits");
    hist->GetYaxis()->SetTitle("Normalized");
    }


    // Draw legend
    legend->Draw();

    // Save the plot as an image
    c->SaveAs("num_rechits_combined.png");

    // Clean up
    for (auto* hist : histograms) {
        delete hist;
    }
    delete legend;
    delete c;
}