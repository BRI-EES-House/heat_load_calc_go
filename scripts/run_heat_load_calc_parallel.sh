#!/bin/bash

# 同時に実行する最大ジョブ数
MAX_JOBS=8

# Loop through regions 1 to 8
for region in {1..8}; do
    for file in data_example1.json jjj.json; do
        # Create a unique output directory for each combination
        output_dir="out/region_${region}_${file%.json}"
        mkdir -p "$output_dir"
        # Execute the command in a subshell and send it to the background
        (
            ./heat_load_calc_go -i example/$file --region=$region -o $output_dir --schedule_saved --weather_saved
        ) &

        # 現在のジョブ数を取得
        job_count=$(jobs -p | wc -l)

        # ジョブ数が最大に達したら待機
        if (( job_count >= MAX_JOBS )); then
            wait -n
        fi
    done
done

# 残っているジョブが完了するのを待つ
wait
