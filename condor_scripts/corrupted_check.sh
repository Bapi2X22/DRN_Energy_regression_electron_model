#!/bin/bash

# Usage: ./check_root_corruption.sh /path/to/root/files

folder="$1"
output_file="corrupted_files.txt"

if [ -z "$folder" ] || [ ! -d "$folder" ]; then
  echo "Usage: $0 /path/to/folder"
  exit 1
fi

# Empty or create the output list
> "$output_file"

echo "Scanning for corrupted ROOT files in: $folder"

# Loop through all .root files in the folder
for file in "$folder"/*.root; do
  if [ -f "$file" ]; then
    # Run ROOT in batch mode and check output
    if root -l -b -q "$file" 2>&1 | grep -q "probably not closed"; then
      echo "Corrupted: $file"
      echo "$file" >> "$output_file"
    fi
  fi
done

echo "Scan complete. Corrupted files (if any) are listed in: $output_file"

