import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2500)

#options = VarParsing.VarParsing('standard')
options = VarParsing.VarParsing('analysis')

#options.register('inputFile',
#        "~/",
#        VarParsing.VarParsing.multiplicity.singleton,
#        VarParsing.VarParsing.varType.string,
#        "File containing a list of the EXACT location of the output file  (default = ~/)"
#        )


options.parseArguments()
#options.inputFile = 'root://eoscms//' + options.inputFile
#print(options.inputFile)
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                #    'root://cms-xrd-global.cern.ch//store/data/Run2018A/EGamma/AOD/12Nov2019_UL2018-v2/710000/B43485A1-DA02-2747-8BD5-C1E17313CC27.root'
                                #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/DoublePhoton_Pt-5To300-gun/AODSIM/FlatPU0to70EdalIdealGT_EdalIdealGT_106X_upgrade2018_realistic_v11_L1v1_EcalIdealIC-v2/280000/DA55E6E8-8F95-CE4C-9FB9-8393ECE23A09.root'
                                #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/GluGluHToGG_M-125_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/40000/9DAAC020-4F79-6641-B3C4-A48EBECECED5.root'
                                'root://cms-xrd-global.cern.ch:1094//store/mc/RunIISummer20UL18RECO/GluGluHToGG_M-125_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/40000/418EDE70-2DF3-714A-8436-71F580A9ED86.root'
                                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.nTuplelize = cms.EDAnalyzer('Photon_RefinedRecHit_NTuplizer',
        isMC = cms.bool(True),
        miniAODRun = cms.bool(False),
        useOuterHits = cms.bool(False),
        rhoFastJet = cms.InputTag("fixedGridRhoFastjetAll"),
        photons = cms.InputTag("gedPhotons"),
        genParticles = cms.InputTag("genParticles"),
        #MVA Based Id
        eleMediumIdMap = cms.InputTag(""),
        eleTightIdMap = cms.InputTag(""),
        #Calo clusters
        eeebClusters = cms.InputTag("particleFlowEGamma:EBEEClusters:RECO"),
        esClusters = cms.InputTag("particleFlowEGamma:ESClusters:RECO"),

        ebNeighbourXtalMap = cms.FileInPath("EM_GNN_ID/EM_Skimmer/data/EB_xtal_dR0p3_map.root"),
        eeNeighbourXtalMap = cms.FileInPath("EM_GNN_ID/EM_Skimmer/data/EE_xtal_dR0p3_map.root")
	)


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("GammaRecHits_ntuple.root"),
#     fileName = cms.string("Tree_Gamma_ABCD.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.egmPhotonIDSequence*process.nTuplelize)

