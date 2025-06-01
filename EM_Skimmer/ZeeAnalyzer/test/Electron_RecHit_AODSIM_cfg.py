import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)

options = VarParsing.VarParsing('standard')

options.register('inputFile',
        "~/",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File containing a list of the EXACT location of the output file  (default = ~/)"
        )


options.parseArguments()
options.inputFile = 'root://eoscms//' + options.inputFile
print(options.inputFile)
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/270000/D35E3B33-20E4-0145-B6F8-239DAC9AEE66.root'
           # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL18RECO/DoubleElectron_FlatPt-1To300/AODSIM/FlatPU0to70RAW_106X_upgrade2018_realistic_v11_L1v1_ext2-v2/40000/FF912577-7001-DE41-BDB0-6FF06EC01DA6.root'
          'file:/eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/test/miniAOD_new_.root'
#                options.inputFile
                ),
                inputCommands=cms.untracked.vstring(
               # 'drop recoTrackExtrasedmAssociation_muonReducedTrackExtras_*_*'
               'keep *'
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)


# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']


for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("gedGsfElectrons"),
        gsfElectronsDRN = cms.InputTag("gsfElectronsDRN"),
        genParticles = cms.InputTag("genParticles"),
        refinedCluster = cms.bool(False),
        isMC = cms.bool(True), 
        miniAODRun = cms.bool(False),
        #MVA Based Id
	eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight")
        )


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("ElectronRecHits_ntuple_4.root"),
     closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.egmGsfElectronIDSequence*process.nTuplelize)

