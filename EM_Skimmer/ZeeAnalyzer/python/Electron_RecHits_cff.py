import FWCore.ParameterSet.Config as cms

def setupElectronNtuplizer(process):
    process.load("EM_GNN_ID.EM_Skimmer.Electron_RecHits_cfi")
