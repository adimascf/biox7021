#!/usr/bin/env bash
set -euo pipefail

# diret the stderr to the log file
exec 2> "${snakemake_log[0]}"

reads="${snakemake_input[reads]}"
kit1="${snakemake_params[kit1]}"
kit2="${snakemake_params[kit2]}"
outdir1=$(mktemp -d)
outdir2=$(mktemp -d)
output="${snakemake_output[reads]}"

dorado demux --kit-name "$kit1" --emit-fastq --output-dir "$outdir1" "$reads"
dorado demux --kit-name "$kit2" --emit-fastq --output-dir "$outdir2" "$reads"

cat "$outdir1"/*/*/*/*/*.fastq > "$output"
cat "$outdir2"/*/*/*/*/*.fastq >> "$output"

