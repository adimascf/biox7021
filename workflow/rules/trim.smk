REPEAT = int(config.get("repeat", 1))

rule trim_fastplong:
	input:
		reads=get_original_fastqs
	log:
		LOGS / "QC/trimming/fastplong/{sample}.log"
	threads: 8
	resources:
		mem_mb=64000,
		runtime=f"{4 * REPEAT}m",
	conda: 
		ENVS / "fastplong.yaml"
	params:
		nofilter="--disable_quality_filtering --disable_length_filtering"
	output:
		reads=RESULTS / "QC/trimming/fastplong/{sample}.fastq",
		json=RESULTS / "QC/trimming/fastplong/{sample}.json",
		html=RESULTS / "QC/trimming/fastplong/{sample}.html"
	benchmark:
		repeat(BENCHMARK / "QC/trimming/fastplong/{sample}.tsv", REPEAT)
	shell:
		"""
		fastplong -i {input.reads} -o {output.reads} {params.nofilter} --thread {threads} \
				--json {output.json} --html {output.html} 2> {log}
		"""

rule trim_porcechop_abi:
	input:
		reads=get_original_fastqs
	log:
		LOGS / "QC/trimming/porechop_abi/{sample}.log"
	threads: 32
	resources:
		mem_mb=256000,
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	params:
		nochimera="--discard_middle"
	output:
		reads=RESULTS / "QC/trimming/porechop_abi/{sample}.fastq"
	benchmark:
		repeat(BENCHMARK / "QC/trimming/porechop_abi/{sample}.tsv", REPEAT)
	shell:
		"""
		porechop_abi -abi -t {threads} {params.nochimera} -i {input.reads} -o {output.reads} 2> {log}
		"""

rule trim_dorado:
	input:
		reads=get_original_fastqs
	log:
		LOGS / "QC/trimming/dorado/{sample}.log"
	resources:
		mem_mb=32000,
		runtime=f"{15 * REPEAT}m"
	container:
		"docker://nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
	params:
		kit1=lambda wildcards: get_sequencing_kits(wildcards)[0],
		kit2=lambda wildcards: get_sequencing_kits(wildcards)[1],
		output_fq="--emit-fastq"
	output:
		reads=RESULTS / "QC/trimming/dorado/{sample}.fastq"
	benchmark:
		repeat(BENCHMARK / "QC/trimming/dorado/{sample}.tsv", REPEAT)
	shell:
		"""
		tmp_fq=$(mktemp -u).fastq
		(dorado trim {params.output_fq} --sequencing-kit {params.kit1} {input.reads} > $tmp_fq) 2> {log}
		(dorado trim {params.output_fq} --sequencing-kit {params.kit2} $tmp_fq > {output.reads}) 2>> {log} 
		"""

rule barbell_trim:
	input:
		reads=get_original_fastqs
	log:
		LOGS / "QC/trimming/barbell/{sample}.log"
	resources:
		mem_mb=32000,
		runtime=f"{20 * REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	params:
		kit1=lambda wildcards: get_sequencing_kits(wildcards)[0],
		kit2=lambda wildcards: get_sequencing_kits(wildcards)[1]
	output:
		reads= RESULTS / "QC/trimming/barbell/{sample}.fastq"
	benchmark:
		repeat(BENCHMARK / "QC/trimming/barbell/{sample}.tsv", REPEAT)
	script:
		"../scripts/trimming/barbell_kit.sh"
