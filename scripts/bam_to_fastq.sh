#!/usr/bin/env bash
set -euo pipefail

files=$1
outdir=$2

for file in `ls $files/*bam`; do
	base=$(basename $file)
	name=$(echo $base | awk -v FS="." '{print $1}')
	echo "converting $name..."
	samtools fastq $file | gzip > $outdir/$name.fastq.gz
done
