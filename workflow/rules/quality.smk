REPEAT = int(config.get("repeat", 1))

tool = "fastplong_quality"
rule quality_fastplong:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / f"QC/quality/{tool}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{5 * REPEAT}m"
	conda:
		ENVS / "fastplong.yaml"
	params:
		no_adapt_trimming="--disable_adapter_trimming"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.fastq"),
		json=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.json"),
		html=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.html")
	benchmark:   
		repeat(BENCHMARK / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.tsv", REPEAT)
	shell:    
		"""   
		fastplong -i {input.reads} -o {output.reads} {params.no_adapt_trimming} \
				--json {output.json} --html {output.html} --verbose 2> {log} 
		"""    
    
tool = "filtlong"    
rule quality_filtlong:    
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / f"QC/quality/{tool}/{{model}}/{{sample}}.log" 
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{60 * REPEAT}m"
	conda:
		ENVS / "filtlong.yaml"
	params:
		min_len="--min_length 1000",
		min_q="--min_mean_q 10",
		pct_keep="--keep_percent 90"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.tsv", REPEAT)
	shell:
		"""
		filtlong {params.min_len} {params.min_q} {params.pct_keep} {input.reads} > {output.reads} 2> {log}
		"""

tool = "chopper_standard"
rule quality_chopper_standard:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / f"QC/quality/{tool}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		standard_filter="-q 20 -l 1000"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.tsv", REPEAT)
	shell:
		"""
		chopper {params.standard_filter} -i {input.reads} > {output.reads} 2> {log}
		"""

tool = "chopper_trim"
rule quality_chopper_trim:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / f"QC/quality/{tool}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		approach="--trim-approach trim-by-quality --cutoff 20 -l 1000"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.tsv", REPEAT)
	shell:
		"""
		chopper {params.approach} -i {input.reads} > {output.reads} 2> {log}
		"""

tool = "chopper_split"
rule quality_chopper_split:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / f"QC/quality/{tool}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		approach="--trim-approach split-by-low-quality --cutoff 20 -l 1000"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.tsv", REPEAT)
	shell:
		"""
		chopper {params.approach} -i {input.reads} > {output.reads} 2> {log}
		"""

tool = "untrimmed"
rule trim_notrim:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
	log:
		LOGS / f"QC/quality/{tool}/{{model}}/{{sample}}.log"
	resources:
		mem="8GiB",
		runtime="10m"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}/{{model}}/{{sample}}.{tool}.fastq")
	shell:
		"""
		# here we just copy the input file, no quality trimming/filtering process being perfomed
		cp {input.reads} {output.reads}
		"""
