#!/bin/bash

for i in {1..5000}; do
    # Define the directory and filename
    directory="Data_1"
    file="${directory}/D1step4_${i}.root"
    output_file="missing_files2.txt"
    search_file="/eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/Dataset2/file_2_copy.txt"
    
    # Check if the file exists
    if [ -f "$file" ]; then
        sed -i "/step4_${i}.root/d" "$search_file"
        continue  # Skip to the next iteration
    fi
done
