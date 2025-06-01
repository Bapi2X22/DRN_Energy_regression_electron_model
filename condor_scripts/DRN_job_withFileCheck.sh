#!/bin/bash

if [ -z $1 ] ; then
  echo "Please use: ./runDRN_final.sh [foldername] [datasetname] [start_idx] [idxcount]" && exit 1;
fi
if [ -z $2 ] ; then
  echo "Please use: ./runDRN_final.sh [foldername] [datasetname] [start_idx] [idxcount]" && exit 1;
fi
if [ -z $3 ] ; then
  echo "Please use: ./runDRN_final.sh [foldername] [datasetname] [start_idx] [idxcount]" && exit 1;
fi
if [ -z $4 ] ; then
  echo "Please use: ./runDRN_final.sh [foldername] [datasetname] [start_idx] [idxcount]" && exit 1;
fi
#mkdir -p $PWD/triton_models
#cp -r /eos/user/s/sosaha/EGM_test_Sonic/CMSSW_13_3_3/src/RecoEgamma/EgammaPhotonProducers/data/models/photonObjectCombined $PWD/triton_models
cp -r /eos/user/b/bbapi/Energy_regression_by_electron_model/CMSSW_13_3_3/src/test/RecoEgamma-EgammaPhotonProducers/models/ /tmp
singularity exec  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:22.07-py3-geometric \
tritonserver --model-repository /tmp/models/ --http-port 8000 --grpc-port 8001 --metrics-port 8002 --allow-http=1 > triton.log 2>&1 & 
#tritonserver --model-repository $PWD/triton_models --http-port 8000 --grpc-port 8001 --metrics-port 8002 --allow-http=1 > triton.log 2>&1 & 
cd /eos/user/b/bbapi/Energy_regression_by_electron_model/CMSSW_13_3_3/src/condor_scripts/
export HOME=/afs/cern.ch/user/b/bbapi
source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsenv;

#====================================== Both noise and energy threshold varied ===================================
#dataset="pfTL235_pedTL235"

#"dataset=pfThresUL18_pedTL235"
IFS=$'\n' files=($(dasgoclient --query="file dataset=$2 instance=prod/phys03"))
unset IFS

if [ ! -f DRN_reg_final_cfg.py ]; then
    echo "Error: DRN_reg_final_cfg.py not found!" >&2
    ls -l  # List directory contents for debugging
    exit 1
fi
for i in "${files[@]:$3:$4}"
do
  echo "Processing file: ${i}"
  # Uncomment and modify the line below if needed
  # echo "Running Skimmer on $folder/${i}_AToGG_RECO_M1000.0.root"
  cmsRun DRN_reg_final_cfg.py inputFile="${i}" datasetname="$1"
  echo "Skimming done"
done
#==============================================================================================
#echo "Checking for files smaller than 200KB in the output folder..."
#output_folder="/eos/user/s/sosaha/EGM_test_Sonic/CMSSW_13_3_3/src/Skimmed"
#small_files=()
#corrupted_files=()
#for i in "${files[@]:$3:$4}"
#do
#  # Extract the file name from the input file path
#  base_name=$(basename "$i")  # This will give something like step4_1171.root
#  output_file="${output_folder}/$1/${base_name}"
  
#  if [ -f "$output_file" ]; then  # Check if the file exists
#    file_size=$(stat -c%s "$output_file")  # Get file size in bytes
#    if (( file_size < 200 * 1024 )); then  # 200KB = 200 * 1024 bytes
#      small_files+=("$i")  # Store the original input file path for reprocessing
#    fi
#    if root -l -b -q "$output_file" 2>&1 | grep -q "probably not closed"; then
#      corrupted_files+=("$i")
#    fi
#  else
#    echo "Warning: Output file $output_file not found. Skipping size check for this file."
#  fi
#done

## If there are small files, process them again
#if [ ${#small_files[@]} -gt 0 ]; then
#  echo "Found ${#small_files[@]} files in the output folder smaller than 200KB. Reprocessing..."
#  for i in "${small_files[@]}"
#  do
#    echo "Reprocessing file: ${i}"
#    cmsRun DRN_reg_final_cfg.py inputFile="${i}" datasetname="$1"
#    echo "Reprocessing done"
#  done
#else
#  echo "No files are smaller than 200KB in the output folder. All done!"
#fi

#if [ ${#corrupted_files[@]} -gt 0 ]; then
#  echo "Found ${corrupted_files[@]} files in the output folder improperly closed. Reprocessing..."
#  for i in "${corrupted_files[@]}"
#  do  
#    echo "Reprocessing file: ${i}"
#    cmsRun DRN_reg_final_cfg.py inputFile="${i}" datasetname="$1"
#    echo "Reprocessing done"
#  done
#else
#  echo "No corrupted files. All done!"
#fi

