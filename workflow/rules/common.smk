from pathlib import Path

def get_original_fastqs(wildcards):
	reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
	return (reads_dir / f"{wildcards.sample}.fastq")

def get_reference_genome(wildcards):
	ref_dir = Path(pep.get_sample(wildcards.sample)["reference_path"])
	return (ref_dir / f"{wildcards.sample}.fna")

def get_mutreference_genome(wildcards):
	ref_dir = Path(pep.get_sample(wildcards.sample)["mutreference_path"])
	return (ref_dir / f"{wildcards.sample}.mutref.fna")

def get_truth_vcf(wildcards):
	ref_dir = Path(pep.get_sample(wildcards.sample)["truth_vcf_path"])
	return (ref_dir / f"{wildcards.sample}.truth.vcf.gz")

def get_sequencing_kits(wildcards):
	kits =  pep.get_sample(wildcards.sample)["sequencing_kits"].split(";")
	return kits[0], kits[1]

