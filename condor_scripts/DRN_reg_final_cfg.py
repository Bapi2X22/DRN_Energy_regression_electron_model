import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import socket
process = cms.Process("DRNregression")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'130X_mcRun3_2023_realistic_v14','')
process.GlobalTag = GlobalTag(process.GlobalTag,'auto:run2_mc','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(24)
                         
options = VarParsing.VarParsing('analysis')
options.register('inputFile',
        "~/",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File containing a list of the EXACT location of the output file  (default = ~/)"
        )
options.register('datasetname',"~/", VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "Folder with name of dataset to store output file  (default = ~/)"
        )
options.parseArguments()
datasetName = str(options.datasetname).split('/')[-1]
infilename = str(options.inputFile).split('/')[-1]
options.inputFile = 'root://cms-xrd-global.cern.ch//' + options.inputFile
print(options.inputFile)                      
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/user/bmarzocc/ECAL_GNN_Regression/FourElectronsGunPt1-500_13TeV-pythia8_RunIISummer20UL18_pfThresTL235_pedTL235_AODSIM/AODSIM/240522_123019/0004/step4_4739.root'
        options.inputFile
    ),
    secondaryFileNames = cms.untracked.vstring()
) 

########### Activate Run 3 2022 IDs [Might need change to the 2023 recommendation, but none exists so far ########
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#dataFormat = DataFormat.MiniAOD ## DataFormat.AOD while running on AOD
dataFormat = DataFormat.AOD ## DataFormat.AOD while running on AOD
switchOnVIDElectronIdProducer(process, dataFormat)
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff']
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.output = cms.OutputModule("PoolOutputModule",
                                   splitLevel = cms.untracked.int32(0),
                                   outputCommands = cms.untracked.vstring("keep *"),
                                   #fileName = cms.untracked.string("DRN_DYtoLL_file1.root")
                                   fileName = cms.untracked.string(options.datasetname+'/'+infilename)
)


from Geometry.CaloEventSetup.CaloGeometryBuilder_cfi import *
CaloGeometryBuilder.SelectedCalos = ['HCAL', 'ZDC', 'EcalBarrel', 'EcalEndcap', 'EcalPreshower', 'TOWER'] # Why is this needed?

#==================================== Running DRN ====================================

process.load("HeterogeneousCore.SonicTriton.TritonService_cff")
process.load("Configuration.ProcessModifiers.photonDRN_cff")
process.load("RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff")
#process.load("PhysicsTools.PatAlgos.slimming.slimming_cff")
#process.load("PhysicsTools.PatAlgos.slimming.patElectronDRNCorrector_cfi")
#process.load("PhysicsTools.PatAlgos.slimming.patPhotonDRNCorrector_cfi")

from Configuration.ProcessModifiers.photonDRN_cff import _photonDRN
#from PhysicsTools.PatAlgos.slimming.patPhotonDRNCorrector_cfi import patPhotonsDRN
from PhysicsTools.PatAlgos.slimming.gsfElectronDRNCorrector_cfi import gsfElectronsDRN

lxplusnode=socket.gethostname()
process.TritonService.servers.append(
    cms.PSet(
        name = cms.untracked.string("local_triton"),
        address = cms.untracked.string(lxplusnode),
        port = cms.untracked.uint32(8001),
        useSsl = cms.untracked.bool(False),
        rootCertificates = cms.untracked.string(""),
        privateKey = cms.untracked.string(""),
        certificateChain = cms.untracked.string(""),)
)
#==================================== Running Ntuplizer =====================================
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff")
#process.load("RecoEgamma.EgammaElectronProducers.gsfElectronDRNCorrector_cfi")  # Explicitly load gedGsfElectrons
#process.load("RecoEgamma.EgammaElectronProducers.gedGsfElectronSequence_cff")
process.load("PhysicsTools.PatAlgos.slimming.gsfElectronDRNCorrector_cfi")

process.nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("gedGsfElectrons"),
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
     #fileName = cms.string("hist.root"),
     #fileName = cms.string("ElectronRecHits_ntuple.root"),
     fileName = cms.string(options.datasetname+'/'+infilename),
     closeFileFast = cms.untracked.bool(True)
  )

#process.p=cms.Path(process.patPhotonsDRN)
#process.p=cms.Path(process.gsfElectronsDRN+process.egmGsfElectronIDSequence*process.nTuplelize)
process.p=cms.Path(process.gsfElectronsDRN)
process.e = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p,process.e)
#process.schedule = cms.Schedule(process.p,process.output)

