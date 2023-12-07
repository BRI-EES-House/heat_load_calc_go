#!/bin/bash

filename1="result_detail_a.csv"
filename2="region_6_jjj/result_detail_a.csv"
echo "Check $filename1 and $filename2"
python3 diff.py "$filename1" "$filename2" --error 0.000001

if [ $? -ne 0 ]; then
    echo "Error occurred while running diff.py. Exiting."
    exit 1
fi
