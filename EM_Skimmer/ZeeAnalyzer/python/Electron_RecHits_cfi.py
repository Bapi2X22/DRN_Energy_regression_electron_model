import FWCore.ParameterSet.Config as cms

nTuplize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        isMC = cms.bool(False),
        miniAODRun = cms.bool(False),
        refinedCluster = cms.bool(True), # use refined clusters by default 
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("gedGsfElectrons"),
        genParticles = cms.InputTag(""),
        #MVA Based Id
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID_Fall17_iso_V2_wp90"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID_Fall17_iso_V2_wp80")
	)
