import FWCore.ParameterSet.Config as cms

nTuplize = cms.EDAnalyzer('Photon_RefinedRecHit_NTuplizer',
        
        isMC = cms.bool(False),
        miniAODRun = cms.bool(False),
        useOuterHits = cms.bool(False),
        rhoFastJet = cms.InputTag("fixedGridRhoFastjetAll"),
        photons = cms.InputTag("gedPhotons"),
        genParticles = cms.InputTag("genParticles"),
        
        # MVA Based Id
        eleMediumIdMap = cms.InputTag(""),
        eleTightIdMap = cms.InputTag(""),

        # xtal maps
        ebNeighbourXtalMap = cms.string("EM_Skimmer/data/EB_xtal_dR0p3_map.xml"),
        eeNeighbourXtalMap = cms.string("EM_Skimmer/data/EE_xtal_dR0p3_map.xml")
	)
