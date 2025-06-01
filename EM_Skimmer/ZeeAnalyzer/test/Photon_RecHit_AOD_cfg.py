import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2500)

options = VarParsing.VarParsing('analysis')
options.parseArguments()
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL18RECO/DoubleElectron_FlatPt-1To300/AODSIM/FlatPU0to70RAW_106X_upgrade2018_realistic_v11_L1v1_ext2-v2/40000/FF912577-7001-DE41-BDB0-6FF06EC01DA6.root'
           #'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/DoublePhoton_FlatPt-5To300/AODSIM/FlatPU0to70RAW_106X_upgrade2018_realistic_v11_L1v1_ext1-v2/00000/6F7BDF5D-2A6A-7342-993F-E8DBABF920C8.root'
   #         'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/250000/35BAE18A-E8EA-254E-B693-9212D3725212.root'
     #       'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/250000/26B3A9F9-290D-6346-8862-2BE71645C57D.root'
      #          options.inputFiles
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)


# add ntuplizer
from EM_GNN_ID.EM_Skimmer.Photon_RecHits_cff import setupNtuplizer
setupNtuplizer(process)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string("GammaRecHits_ntuple.root"),
     closeFileFast = cms.untracked.bool(True)
  )

process.p = cms.Path(process.egmPhotonIDSequence*process.nTuplize)
