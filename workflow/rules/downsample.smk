rule faidx_self:
	input:
		reference=get_reference_genome,
	log:
		LOGS / "faidx_reference/{sample}.log",
	resources:
		mem="1GiB",
		runtime="5m"
	conda:
		ENVS / "rasusa.yaml"
	output:
		fasta=RESULTS / "reference/{sample}.fa",
		faidx=RESULTS / "reference/{sample}.fa.fai",
	shell:
		"""
		cp {input.reference} {output.fasta} 2> {log}
		samtools faidx {output.fasta} 2>> {log}
		"""

use rule faidx_self as faidx_mutref with: 
	input:
		reference=get_mutreference_genome
	log:
		LOGS / "faidx_reference/{sample}.mutref.log"
	output:
		fasta=RESULTS / "mutreference/{sample}.mutref.fa",
		faidx=RESULTS / "mutreference/{sample}.mutref.fa.fai"

rule downsample_rasusa:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{sample}.{trimmer}.fastq",
		faidx=rules.faidx_self.output.faidx
	log:
		LOGS / "QC/downsampling/{trimmer}/{sample}.{trimmer}.rasusa.log"
	resources:
		mem="64GiB",
		runtime="30m"
	conda:
		ENVS / "rasusa.yaml"
	params:
		depth="100",
		seed="1"
	output:
		reads=temp(RESULTS / "QC/downsampling/{trimmer}/{sample}.{trimmer}.rasusa.fastq"),
	shell:
		"""
		rasusa reads -c {params.depth} -g {input.faidx} -s {params.seed} -o {output.reads} {input.reads} 2> {log}
		"""
