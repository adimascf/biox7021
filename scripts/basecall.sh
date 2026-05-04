#!/usr/bin/env bash
set -euo pipefail


dir=$1
model=$2

for run in "$dir"/*/; do
    name=$(basename "${run%/}")
	dorado basecaller --no-trim --recursive $model $run > $dir/$name.bam
done

# ssubmit -t 1d -m 16g basecall "./scripts/basecall.sh <inputdir> <model>" -- --qos gpu -p gpu_cuda --gres=gpu:h100:1
