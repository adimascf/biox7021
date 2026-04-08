REPEAT = int(config.get("repeat", 1))

trimmer = "fastplong"
rule trim_fastplong:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{4 * REPEAT}m",
	conda:
		ENVS / "fastplong.yaml"
	params:
		nofilter="--disable_quality_filtering --disable_length_filtering"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq"),
		json=RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.json",
		html=RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.html"
	benchmark:
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell:
		"""
		fastplong -i {input.reads} -o {output.reads} {params.nofilter} --thread {threads} \
				-s TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT -e ACGTAACTGAACGAAGTACAGG \
				--json {output.json} --html {output.html} --verbose 2> {log}
		"""

trimmer = "porechop_abi_split"
# by default, it will split reads when an adapter is found in the middle
rule trim_porcechop_abi_split:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	threads: 32 
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	benchmark: 
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell: 
		"""
		porechop_abi -abi -t {threads} -i {input.reads} -o {output.reads} 2> {log}
		"""

trimmer = "porechop_abi_discard"
rule trim_porcechop_abi_discard:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	threads: 32
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	params:
		nochimera="--discard_middle"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	benchmark: 
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell: 
		"""
		porechop_abi -abi -t {threads} {params.nochimera} -i {input.reads} -o {output.reads} 2> {log}
		"""

trimmer = "porechop_abi_nocheck"
rule trim_porcechop_abi_nocheck:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	threads: 32
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	params:
		nocheck="--no_split" # will skip chimera searching
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	benchmark: 
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell:  
		""" 
		porechop_abi -abi -t {threads} {params.nocheck} -i {input.reads} -o {output.reads} 2> {log}
		"""  

trimmer = "dorado"
rule trim_dorado:  
	input:  
		reads=get_original_fastqs  
	log:  
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"  
	resources:  
		mem="32GiB",  
		runtime=f"{15 * REPEAT}m"  
	container:  
		"docker://nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
	params:  
		kit1=lambda wildcards: get_sequencing_kits(wildcards)[0],  
		kit2=lambda wildcards: get_sequencing_kits(wildcards)[1],  
		output_fq="--emit-fastq"  
	output:  
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	benchmark:  
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell: 
		"""
		tmp_fq=$(mktemp -u).fastq
		(dorado trim {params.output_fq} --sequencing-kit {params.kit1} {input.reads} > $tmp_fq) 2> {log}
		(dorado trim {params.output_fq} --sequencing-kit {params.kit2} $tmp_fq > {output.reads}) 2>> {log} 
		"""

trimmer = "dorado_demux"
rule trim_dorado_demux:  
	input:  
		reads=get_original_fastqs  
	log:  
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"  
	resources:  
		mem="32GiB",  
		runtime=f"{12 * REPEAT}h"  
	container:  
		"docker://nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
	params:  
		kit1=lambda wildcards: get_sequencing_kits(wildcards)[0],  
		kit2=lambda wildcards: get_sequencing_kits(wildcards)[1],  
	output:  
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}_raw.fastq")
	benchmark:  
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	script: 
		"../scripts/trimming/dorado_demux.sh"

rule seqkit_dedup:
	input:
		reads=rules.trim_dorado_demux.output.reads
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}_seqkit.log"
	resources:
		mem="32GiB",
		runtime=f"{30*REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	output:
		reads=RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq" 
	shell:
		"""
		# remove duplicates by read ID
		seqkit rmdup -o {output.reads} {input.reads} 2> {log}
		"""

trimmer = "barbell_default"
rule barbell_trim:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	resources:
		mem="32GiB",
		runtime=f"{20 * REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	params:
		kit1=lambda wildcards: get_sequencing_kits(wildcards)[0],
		kit2=lambda wildcards: get_sequencing_kits(wildcards)[1],
		maximize="no"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	script:
		"../scripts/trimming/barbell_kit.sh"

trimmer = "barbell_max"
rule barbell_trim_max:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	resources:
		mem="32GiB",
		runtime=f"{20 * REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	params:
		kit1=lambda wildcards: get_sequencing_kits(wildcards)[0],
		kit2=lambda wildcards: get_sequencing_kits(wildcards)[1],
		maximize="yes"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.tsv", REPEAT)
	script:
		"../scripts/trimming/barbell_kit.sh"

trimmer = "untrimmed"
rule trim_notrim:
	input:
		reads=get_original_fastqs
	log:
		LOGS / f"QC/trimming/{trimmer}/{{sample}}.log"
	resources:
		mem="8GiB",
		runtime="10m"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{sample}}.{trimmer}.fastq")
	shell:
		"""
		# here we just copy the input file, no trimming process being perfomed
		cp {input.reads} {output.reads}
		"""
