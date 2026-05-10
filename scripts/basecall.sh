#!/usr/bin/env bash
set -euo pipefail


dir=$1
model=$2
outdir=$3
trim_flag=${4:-""} # Default is empty (trimming enabled)

mkdir -p "$outdir"

for run in "$dir"/*/; do
    name=$(basename "${run%/}")
	dorado basecaller $trim_flag --recursive "$model" "$run" > "$outdir/$name.bam"
done

# ssubmit -t 1d -m 16g basecall "./scripts/basecall.sh <inputdir> <model>" -- --qos gpu -p gpu_cuda --gres=gpu:h100:1
