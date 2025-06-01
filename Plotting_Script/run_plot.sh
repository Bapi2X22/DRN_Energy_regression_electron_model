#!/bin/bash

# Loop over runlist indices
for i in {1..4}; do
    # Generate the runlist file name
    runlist="runlist_${i}.txt"
    
    # Loop over positions for '1' in the array
    for j in {0..4}; do
        # Generate the binary array
        params=(0 0 0 0 0)
        params[$j]=1
        
        # Determine etaX value
        eta_index=$((j + 1))
        
        # Generate the output file name
        output_file="Plot_R9_eta${eta_index}_${i}.root"
        
        # Run the command
        ./analyzeHGCMuons_D1 "$runlist" "$output_file" data 2500 "${params[@]}"
    done
done