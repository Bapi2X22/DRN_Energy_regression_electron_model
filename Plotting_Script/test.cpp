#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <vector>

using namespace std;

void printBranch(const char* filename, const char* treeName, const char* branchName) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }

    // Get the TTree
    TTree *tree = (TTree*)file->Get(treeName);
    if (!tree) {
        cerr << "Error: Tree " << treeName << " not found in file " << filename << endl;
        return;
    }

    cout << "Tree '" << treeName << "' found in file '" << filename << "'\n";
    cout << "Branch '" << branchName << "' entries:\n";

    // Check if the branch exists
    TBranch *branch = tree->GetBranch(branchName);
    if (!branch) {
        cerr << "Error: Branch " << branchName << " not found in tree " << treeName << endl;
        return;
    }

    // Attempt to handle both scalar and vector types
    vector<float> *vecValue = nullptr;
    float scalarValue;

    // Check if the branch is a vector (dynamic)
    tree->SetBranchAddress(branchName, &vecValue);
    if (tree->GetEntry(0) && vecValue) {
        cout << "[Detected vector<float> branch]\n";
        for (Long64_t i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);
            cout << "Entry " << i << ": ";
            for (size_t j = 0; j < vecValue->size(); j++) {
                cout << vecValue->at(j) << " ";
            }
            cout << endl;
        }
    } else {
        // If not a vector, assume it's a scalar (e.g., float)
        tree->SetBranchAddress(branchName, &scalarValue);
        cout << "[Detected scalar branch]\n";
        for (Long64_t i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);
            cout << "Entry " << i << ": " << scalarValue << endl;
        }
    }

    // Cleanup
    file->Close();
}
