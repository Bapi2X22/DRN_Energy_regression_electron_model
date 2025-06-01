import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import socket
import subprocess
import os
import time
options = VarParsing.VarParsing('analysis')

process = cms.Process("ZeeDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

# from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag,'auto:run2_mc','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10)

options.register('inputFile',
                '',  # Default value (empty)
                 VarParsing.VarParsing.multiplicity.singleton,  # Single input
                 VarParsing.VarParsing.varType.string,  # String type
                 "Input file for the cmsRun job")

options.register('datasetname',"~/", VarParsing.VarParsing.multiplicity.singleton,
         VarParsing.VarParsing.varType.string,
         "Folder with name of dataset to store output file  (default = ~/)"
         )

#options.register('outputFile',
#                 'miniAOD_new.root',
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.string,
#                 "Output file name")
# Parse command-line arguments
options.parseArguments()
datasetName = str(options.datasetname).split('/')[-1]
infilename = str(options.inputFile).split('/')[-1]
options.inputFile = 'root://cms-xrd-global.cern.ch//' + options.inputFile
print(options.inputFile)


# Validate the input
#if not options.inputFile:
#    raise ValueError("No input file specified. Use '--inputFile' argument.")


process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
       # 'root://cms-xrd-global.cern.ch//store/mc//RunIISummer20UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/130000/9E0A8A10-D7B9-3C48-A765-9D54EF930E38.root'
      # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/270004/22C4E1FC-D896-2345-B4F5-05CAD89E343E.root'
     # 'root://cms-xrd-global.cern.ch//store/user/bmarzocc/ECAL_GNN_Regression/FourElectronsGunPt1-500_13TeV-pythia8_RunIISummer20UL18_pfThresTL235_pedTL235_AODSIM/AODSIM/240522_123019/0000/step4_513.root'
      options.inputFile
#      'root://cms-xrd-global.cern.ch//store/user/bmarzocc/ECAL_GNN_Regression/FourElectronsGunPt1-500_13TeV-pythia8_RunIISummer20UL18_pfThresTL235_pedTL235_AODSIM/AODSIM/240522_123019/0000/step4_680.root'
#       'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/DoublePhoton_Pt-5To300-gun/AODSIM/FlatPU0to70EdalIdealGT_EdalIdealGT_106X_upgrade2018_realistic_v11_L1v1_EcalIdealIC-v2/130000/CA62D711-F4ED-3741-AA56-5D344BD955AF.root'
#       'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/DoublePhoton_Pt-5To300-gun/AODSIM/FlatPU0to70EdalIdealGT_EdalIdealGT_106X_upgrade2018_realistic_v11_L1v1_EcalIdealIC-v2/130000/EFE5A23C-7B5C-CF48-B803-97E404663C12.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
) 

#file_index = os.environ.get('OUT_INDEX', '1')

######################Activate Run 3 2022 IDs [Might need change to the 2023 recommendation, but none exists so far]##########################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
#my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#for idmod in my_id_modules:
#        setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)


# process.output = cms.OutputModule("PoolOutputModule",
#                                    splitLevel = cms.untracked.int32(0),
#                                    outputCommands = cms.untracked.vstring("keep *"),
#                                    #fileName = cms.untracked.string("DRN_DYtoLL_file1.root")
#                                    fileName = cms.untracked.string(options.datasetname+'/'+infilename)
# )
##############################################################################################################################################################

########################## Make Photon regressed energies and the IDs accessible from the electron pointer ########################################### 

#################################################################################################################################

#process.load('ScaleAndSmearingTools.Dumper.Zee_dumper_MINIAOD_cfi') # Runs the ele energy producer and sets up the dumper
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string("output.root")
#)

#process.output = cms.OutputModule("PoolOutputModule",
#                                   splitLevel = cms.untracked.int32(0),
#                                   outputCommands = cms.untracked.vstring("keep *"),
#                                   fileName = cms.untracked.string("miniAOD_New.root")
#                          #      fileName = cms.untracked.string(options.outputFile)
#)




from Geometry.CaloEventSetup.CaloGeometryBuilder_cfi import *
CaloGeometryBuilder.SelectedCalos = ['HCAL', 'ZDC', 'EcalBarrel', 'EcalEndcap', 'EcalPreshower', 'TOWER'] # Why is this needed?

#process.eleNewEnergies_step = cms.Path(process.egmGsfElectronIDSequence+process.eleNewEnergiesProducer+process.slimmedECALELFElectrons+process.zeedumper)
process.load("HeterogeneousCore.SonicTriton.TritonService_cff")
#process.load("PhysicsTools.PatAlgos.slimming.patPhotonDRNCorrector_cfi")
#process.load("PhysicsTools.PatAlgos.slimming.gedPhotonDRNCorrector_cfi")
process.load("Configuration.ProcessModifiers.photonDRN_cff")
process.load("RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff")
#process.load("PhysicsTools.PatAlgos.slimming.patElectronDRNCorrector_cfi")
process.load("PhysicsTools.PatAlgos.slimming.gsfElectronDRNCorrector_cfi")

from Configuration.ProcessModifiers.photonDRN_cff import _photonDRN
#from PhysicsTools.PatAlgos.slimming.patPhotonDRNCorrector_cfi import patPhotonsDRN
#from PhysicsTools.PatAlgos.slimming.gedPhotonDRNCorrector_cfi import gedPhotonsDRN
from PhysicsTools.PatAlgos.slimming.gsfElectronDRNCorrector_cfi import gsfElectronsDRN
from RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_tools import *

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')



#    print(triton_process.stderr.read().decode())

process.TritonService.servers.append(
    cms.PSet(
        name = cms.untracked.string("local_triton"),
        address = cms.untracked.string(socket.gethostname()),
        port = cms.untracked.uint32(9001),
        useSsl = cms.untracked.bool(False),
        rootCertificates = cms.untracked.string(""),
        privateKey = cms.untracked.string(""),
        certificateChain = cms.untracked.string(""),)
)
process.TritonService.verbose = cms.untracked.bool(True)

#from Configuration.ProcessModifiers.photonDRN_cff import _photonDRN
#from PhysicsTools.PatAlgos.slimming.patPhotonDRNCorrector_cfi import patPhotonsDRN
#process.dumper_step = cms.Path(process.zeedumper)
#process.out = cms.EndPath(process.output)
#process.p=cms.Path(process.gedPhotonsDRN+gsfElectronsDRN)
input_string = options.inputFile
output = os.path.basename(input_string)

process.nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("gedGsfElectrons"),
#        photons = cms.InputTag("gedPhotons"),
#        gedPhotonsDRN = cms.InputTag("gedPhotonsDRN"),
        gsfElectronsDRN = cms.InputTag("gsfElectronsDRN"),
        genParticles = cms.InputTag("genParticles"),
        refinedCluster = cms.bool(False),
        isMC = cms.bool(True),
        miniAODRun = cms.bool(False),
        #MVA Based Id
#        eleLooseIdMap = cms.InputTag(""),
#        eleMediumIdMap = cms.InputTag(""),
#        eleTightIdMap = cms.InputTag(""),
        eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight")
        )

#process.TFileService = cms.Service("TFileService",
#     #fileName = cms.string("hist.root"),
#     #fileName = cms.string("ElectronRecHits_ntuple.root"),
#     fileName = cms.string(output),
#     closeFileFast = cms.untracked.bool(True)
#  )

process.TFileService = cms.Service("TFileService",
     #fileName = cms.string("hist.root"),
     #fileName = cms.string("ElectronRecHits_ntuple.root"),
     fileName = cms.string(options.datasetname+'/'+infilename),
     closeFileFast = cms.untracked.bool(True)
  )

process.p=cms.Path(process.gsfElectronsDRN+process.egmGsfElectronIDSequence*process.nTuplelize)
# process.e = cms.EndPath(process.output)
#process.schedule = cms.Schedule(process.eleNewEnergies_step)
#process.schedule = cms.Schedule(process.eleNewEnergies_step)
process.schedule = cms.Schedule(process.p)
# process.schedule = cms.Schedule(process.p, process.e)
