#!/usr/bin/env bash
set -euo pipefail


dir=$1

for run in "$dir"/*/; do
    name=$(basename "${run%/}")
	dorado basecaller --no-trim --recursive sup $run > $dir/$name.bam
done

# ssubmit -t 1d -m 16g basecall "./scripts/basecall.sh <inputdir>" -- --qos gpu -p gpu_cuda --gres=gpu:h100:1
