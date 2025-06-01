#!/bin/bash

# Set your CMSSW config file name here
CONFIG_FILE="Electron_RecHit_AODSIM_cfg_mod.py"

# Check if input file list is given
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 filnames.txt"
    exit 1
fi

FILELIST="$1"

# Loop through each line (ROOT file) in the file list
while IFS= read -r rootfile; do
    # Skip empty lines or lines starting with '#'
    [[ -z "$rootfile" || "$rootfile" == \#* ]] && continue

    echo "Running cmsRun on $rootfile"
    cmsRun "$CONFIG_FILE" inputFile="file:$rootfile"

done < "$FILELIST"

