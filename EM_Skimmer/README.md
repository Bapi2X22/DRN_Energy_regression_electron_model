## Setup
To get started, run the following lines on lxplus:

```bash
cmsrel CMSSW_10_6_8
cd CMSSW_10_6_8/src/
cmsenv
mkdir EM_GNN_ID
cd EM_GNN_ID
git clone git@github.com:UMN-CMS/EM_Skimmer.git
cd EM_Skimmer
scram b -j8
```

## Test Run
```bash
cd test
cmsRun Photon_RecHit_AOD_cfg.py
```

Note: Max events are set to 1000. Change to -1 before submitting the crab jobs. The datasets to run on editing the crab config are in Gamma_Jet_Dataset.txt

## CRAB
To generate config file for CRAB job submission, create a file in test/configs (i.e. see example_config.yml). Then run the following lines

```bash
source setpath.sh
./make_crab_cfg.py -c <path-to-config-file> -t <production-tag>
```

The script will produce python configuration files to submit CRAB jobs. Run
```bash
source SubmitCRAB.sh
```

to submit all the CRAB jobs at once. If needed, adjust the number of files per unit and number of units per job as per the recommendations using the dryrun option.
The work area for CRAB jobs will be stored in a directory with the corresponding production tag under crab_projects. Each directory will further contain job directories depending on the number of datasets provided as input. 
```bash
crab submit -c <crab-cfg-file> --dryrun
```
