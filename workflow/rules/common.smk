from pathlib import Path

def get_original_fastqs(wildcards):
	reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
	return (reads_dir / f"{wildcards.sample}.fastq.gz")

def get_sequencing_kit(wildcards):
	return pep.get_sample(wildcards.sample)["sequencing_kit"]
