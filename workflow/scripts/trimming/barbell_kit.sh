#!/usr/bin/env bash
set -euo pipefail

# diret the stderr to the log file
exec 2> "${snakemake_log[0]}"

reads="${snakemake_input[reads]}"
kit1="${snakemake_params[kit1]}"
kit2="${snakemake_params[kit2]}"
outdir=$(mktemp -d)
output="${snakemake_output[reads]}"
max="${snakemake_params[maximize]}"

if [[ "$max" == "yes" ]]; then
	 
	# outputing to the same directory because they have different barcode names, 
	# NB prefix for the native kit, RB or BC or RBK prefix for the rapid kit
	# https://nanoporetech.com/document/chemistry-technical-document#barcode-sequences
	barbell kit --maximize --kit "$kit1" --input "$reads" --output "$outdir"
	barbell kit --maximize --kit "$kit2" --input "$reads" --output "$outdir"

	# concatenate all of the "demultiplexed" output files, and the we remove the duplicate by read name
	cat "$outdir"/*.fastq | seqkit rmdup -n -o "$output"

elif [[ "$max" == "no" ]]; then
 
	# outputing to the same directory because they have different barcode names, 
	# NB prefix for the native kit, RB or BC or RBK prefix for the rapid kit
	# https://nanoporetech.com/document/chemistry-technical-document#barcode-sequences
	barbell kit --kit "$kit1" --input "$reads" --output "$outdir"
	barbell kit --kit "$kit2" --input "$reads" --output "$outdir"

	# concatenate all of the "demultiplexed" output files, and the we remove the duplicate by read name
	cat "$outdir"/*.fastq | seqkit rmdup -n -o "$output"
fi

