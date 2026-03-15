rule trim_fastplong:
	input:
		reads=get_original_fastqs
	log:
		LOGS / "QC/trimming/fastplong/{sample}.log"
	threads: 8
	conda: 
		ENVS / "fastplong.yaml"
	params:
		nofilter="--disable_quality_filtering --disable_length_filtering"
	output:
		reads=RESULTS / "QC/trimming/fastplong/{sample}.fastq.gz",
		json=RESULTS / "QC/trimming/fastplong/{sample}.json",
		html=RESULTS / "QC/trimming/fastplong/{sample}.html"
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
	conda:
		ENVS / "porechop_abi.yaml"
	output:
		reads=RESULTS / "QC/trimming/porechop_abi/{sample}.fastq.gz"
	shell:
		"""
		porechop_abi -abi -t {threads} --no_split -i {input.reads} -o {output.reads} 2> {log}
		"""

rule trim_dorado:
	input:
		reads=get_original_fastqs
	log:
		LOGS / "QC/trimming/dorado/{sample}.log"
	container:
		"docker://nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
	params:
		kit_name=get_sequencing_kit,
		output_fq="--emit-fastq"
	output:
		reads=RESULTS / "QC/trimming/dorado/{sample}.fastq.gz"
	shell:
		"""
		(dorado trim {params.output_fq} --sequencing-kit {params.kit_name} {input.reads} | gzip > {output.reads}) 2> {log} 
		"""
