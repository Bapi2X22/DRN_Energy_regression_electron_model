# -------------------------------------
# Script to generate CRAB jobs
# -------------------------------------

#!/usr/bin/env python

import os, sys

from yaml_utils import *
from collections import OrderedDict

from datetime import datetime
from time_utils import *


class crab_job_manager():
  
    def __init__(self, file_ = ''):

        self.config_file = file_
        self.config = None
        self.mc = False
        self.valid_config = False
        self.production_tag = ''
        self.parameter_list = ['InputDataset', 'OutputFile', 'JSON', 'GlobalTag',\
                 'RequestName','WorkArea', 'TransferOutputs','TransferLogs',\
                 'PsetName','MaxJobRunTime','InputDBS','Splitting',\
                 'UnitsPerJob','FilePath','StorageSite','UseWhitelist',\
                 'UseBlacklist','Whitelist','Blacklist', 'IgnoreGlobalBlacklist']

    def check_config(self):

        for key_ in self.parameter_list:
            if (key_=='GlobalTag' or key_=='JSON') and self.mc: continue
            if key_ not in self.config:
                print 'Specify '+key_+' in the config file!'
                return False

        return True

    def read_config_file(self, file_):

        self.config_file = file_
        self.config = ordered_load(open(file_), Loader=yaml.SafeLoader)

        self.config['RequestName'] = ''
        self.config['WorkArea'] = 'crab_projects/{}'.format(self.production_tag)

        if not self.check_config():
            print 'Loading config failed. Check example_config.yml for help.'
            self.valid_config = False

        else: self.valid_config = True

    def produce_crab_submission_script(self):

        if not self.valid_config:
            print 'Please load config file first!'
            return

        cfg_list = []
        for dataset in self.config['InputDataset']:
            
            dataset_tag = dataset.replace('/', '_')
            self.config['RequestName'] = dataset_tag[1:99]
            pyscript_name = 'CRAB_{}{}.py'.format(self.production_tag, dataset_tag)

            with open('crab_template.txt', 'read') as tmp:
                with open(pyscript_name, 'w') as sub:

                        plain_text = tmp.readlines()
                        new_text = []

                        for line in plain_text:
                            newline = line
                            if 'Template file for CRAB job submission' in line:
                                newline = '# CRAB Job: {}\n'.format(print_creation_timestamp())
                            for key_ in self.parameter_list:
                                if not key_ in self.config: continue
                                par_ = self.config[key_]
                                if '<{}>'.format(key_) in line:
                                    # If the key is InputDataset, use only one dataset at a time
                                    if 'InputDataset' in key_:
                                        newline = 'config.Data.inputDataset\t= "{}"\n'.format(dataset)
                                    # Check for the type of the parameter
                                    # If the parameter is a string replace it with a string
                                    # Else if it's a float or an integer or a list, convert to a string
                                    elif type(par_)==str:
                                        newline = line.replace('<{}>'.format(key_), "'{}'".format(par_))
                                    else:
                                        newline = line.replace('<{}>'.format(key_), str(par_))
                                        if key_=='Blacklist' and self.config['UseBlacklist']==False:
                                            newline = '#'+newline
                                        if key_=='Whitelist' and self.config['UseWhitelist']==False:
                                            newline = '#'+newline
                            new_text.append(newline)


                        for line in new_text:
                            sub.write(line)

                        cfg_list.append(pyscript_name)

        return cfg_list

    def setup_crab(self, file_, tag_, mc_):

        self.mc = mc_
        self.production_tag = tag_
        self.read_config_file(file_)
        _ = self.produce_crab_submission_script()

        return _
