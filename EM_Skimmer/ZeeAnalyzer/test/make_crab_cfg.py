#!/usr/bin/env python
#----------------------------------------------------------------
# A simple script to generate CMSSW cfg files for CRAB submission
#----------------------------------------------------------------

from yaml_utils import *
from crab_manager import crab_job_manager
import argparse

parser = argparse.ArgumentParser(description='Generate CRAB submission scripts')
parser.add_argument('-c', '--config', type=str, default='configs/example_config.yml', help='Path to the config file.')
parser.add_argument('-t', '--tag', type=str, default='TestProduction_20_09_2021',\
                    help='Specify the tag for CRAB production. The directories in the crab area will be names accordingly.')
parser.add_argument('-m', '--mc', action="store_true", help='By default this runs on data. Use this option to run on MC.')
args = parser.parse_args()

manager = crab_job_manager()
cfgs = manager.setup_crab(args.config, args.tag, args.mc)

with open('SubmitCRAB.sh', 'w') as sub_:
    for cfg in cfgs:
        sub_.write('crab submit -c {}\n'.format(cfg))
