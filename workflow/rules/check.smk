rule check_quality_ends:
	input:
		reads= RESULTS / "QC/trimming/{trimmer}/{sample}.{trimmer}.fastq"
	resources:
		mem="8GiB",
		runtime="10m"
	log:
		LOGS / "QC/ends/{trimmer}/{sample}.{trimmer}.log"
	output:
		tsv=RESULTS / "QC/ends/{trimmer}/{sample}.{trimmer}.ends.tsv"
	shell:
		"""
		workflow/bin/qends -n 150 -m median -i {input.reads} > {output.tsv} 2> {log}
		"""

rule plot_quality_ends:
	input:
		qends=list(qends_files)
	resources:
		mem="4GiB",
		runtime="5m"
	log:
		LOGS / "QC/ends/trimming_plot_quality_ends.log"
	conda:
		ENVS / "assess_variant_python.yaml"
	output:
		plot=FIGURES / "QC/ends/trimming_plot_quality_ends.png"
	script:
		"../scripts/plot_quality_ends.py"
	
rule check_read_stats:
	input:
		reads=rules.check_quality_ends.input.reads
	resources:
		mem="8GiB",
		runtime="15m"
	log:
		LOGS / "QC/stats/{trimmer}/{sample}.{trimmer}.log"
	conda:
		ENVS / "barbell_seqkit.yaml"
	output:
		tsv=RESULTS / "QC/stats/{trimmer}/{sample}.{trimmer}.stats.tsv"
	shell:
		"""
		seqkit stats {input.reads} > {output.tsv} 2> {log}
		"""
