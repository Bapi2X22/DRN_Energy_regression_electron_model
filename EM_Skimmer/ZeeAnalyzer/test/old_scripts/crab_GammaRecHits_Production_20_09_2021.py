# -------------------------------------
# CRAB Job: Created on 2021/9/21 15:53:11 UTC
# -------------------------------------

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

config.General.requestName = 'GammaRecHits_Production_20_09_2021'
config.General.workArea = 'crab_projects/GammaRecHits_Production_20_09_2021/'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.outputFiles= 'GammaRecHits_ntuple.root'

config.JobType.psetName = 'test/Photon_RecHit_AOD_cfg.py'

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2750
#config.JobType.sendPythonFolder = True

config.Data.inputDataset  = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM'

config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased' #'LumiBased' / 'FileBased'
config.Data.unitsPerJob   = 30 #30000
config.Data.outLFNDirBase = '/store/user/bjoshi'
config.Data.publication   = False

# GRID
config.Site.storageSite   = 'T2_US_Florida'
#config.Site.Whitelist     = ['T3_CH_CERN_CAF']
#config.Site.Blacklist     = ['T1_US_FNAL', 'T2_UA_KIPT', 'T2_UK_London_Brunel', 'T2_CH_CSCS', 'T2_US_*']
config.Site.ignoreGlobalBlacklist = False

crabSubmit(config)
