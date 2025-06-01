try:
  from CRABClient.UserUtilities import config
except ImportError:
  print
  print 'ERROR: Could not load CRABClient.UserUtilities.  Please source the crab3 setup:'
  print 'source /cvmfs/cms.cern.ch/crab3/crab.sh'
  exit(-1)
from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException
#import FWCore.ParameterSet.Config as cms
def crabSubmit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException, hte:
      print '-----> there was a problem. see below.'
      print hte.headers
      print 'quit here'
      exit(-1)
    exit(0)

config = config()

config.General.requestName = ''
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.outputFiles= [
#    'nTupleMC.root'
 #   'nTupleMC_Pt20to40.root'
    'nTupleMC.root'
#    'nTupleMC_Pt40toInf.root'
]

config.JobType.psetName = 'test/Photon_RecHit_AOD_cfg.py'

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2750
#config.JobType.sendPythonFolder = True

#config.Data.inputDataset  = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM'
#config.Data.inputDataset  = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM'
config.Data.inputDataset  = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM'
#config.Data.inputDataset  = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased' #'LumiBased'
config.Data.unitsPerJob   = 10 #30000
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/rchatter/2018_UL_Photon_MC_v3'
config.Data.outLFNDirBase = '/store/user/aevans/FAIR/'
config.Data.publication   = False

# GRID
config.Site.storageSite   =  'T2_CH_CERN'
#config.Site.whitelist     = ['T3_CH_CERN_CAF']
#config.Site.blacklist     = ['T1_US_FNAL','T2_UA_KIPT','T2_UK_London_Brunel','T2_CH_CSCS','T2_US_*']
#config.Site.ignoreGlobalBlacklist = True

crabSubmit(config)
