from pathlib import Path

def get_original_fastqs(wildcards):
	reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
	return (reads_dir / f"{wildcards.sample}.fastq")

def get_sequencing_kits(wildcards):
	kits =  pep.get_sample(wildcards.sample)["sequencing_kits"].split(";")
	return kits[0], kits[1]
