#!/bin/bash

# Loop through regions 1 to 8
for region in {1..8}; do
    # Execute the command with each region and both json files
    for file in data_example1.json jjj.json; do
        # Create a unique output directory for each combination
        output_dir="out/region_${region}_${file%.json}"
        mkdir -p "$output_dir"
        # Execute the command
        ./heat_load_calc_go -i example/$file --region=$region -o $output_dir --schedule_saved --weather_saved
    done
done