rule make_full_bed:
	input:
		faidx=rules.faidx_mutref.output.faidx
	log:
		LOGS / "make_fule_bed/{sample}.bed.log"
	resources:
		mem_mb=500,
		runtime="1m"
	output:
		bed=RESULTS / "mutreference/{sample}.mutref.bed"
	shell:
		"""
		awk '{{print $1"\t"0"\t"$2}}' {input.faidx} | sort -k1,1 -k2,2n > {output.bed} 2> {log}
		"""

rule assess_calls_vcfdist:
	input:
		query_vcf=RESULTS / "calling/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.clair3.vcf.gz",
		truth_vcf=get_truth_vcf,  
		mutreference=get_mutreference_genome,
		faidx=rules.faidx_mutref.output.faidx,
		bed=rules.make_full_bed.output.bed
	log:
		LOGS / "assess_mutref_calls/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.log"
	resources:
		mem="32GiB",
		runtime="30m"
	params:
		opts="--largest-variant 50 --credit-threshold 1.0",
		prefix=lambda wildcards: RESULTS / f"assess/call/{wildcards.tool}-{wildcards.trimmer}/{wildcards.depth}x/{wildcards.model}/{wildcards.sample}"
	output:  
		pr_summary=RESULTS / "assess/call/{tool}-{trimmer}/{depth}x/{model}/{sample}.precision-recall-summary.tsv",
		pr=RESULTS / "assess/call/{tool}-{trimmer}/{depth}x/{model}/{sample}.precision-recall.tsv",
		summary=RESULTS / "assess/call/{tool}-{trimmer}/{depth}x/{model}/{sample}.summary.vcf",
		query=RESULTS / "assess/call/{tool}-{trimmer}/{depth}x/{model}/{sample}.query.tsv",
		truth=RESULTS / "assess/call/{tool}-{trimmer}/{depth}x/{model}/{sample}.truth.tsv"
	container:     
		"docker://timd1/vcfdist:v2.6.4"
	shell: 
		"""
		exec 2> {log}
		echo "Calculated maximum QUAL score..." 1>&2
		MAX_QUAL=$(bgzip -dc {input.query_vcf} | grep -v '^#' | cut -f 6 | sort -gr | sed -n '1p')
		echo "MAX_QUAL=$MAX_QUAL" 1>&2
		echo "Running vcfdist..." 1>&2
		vcfdist {input.query_vcf} {input.truth_vcf} {input.mutreference} {params.opts} \
				-p {params.prefix}. -b {input.bed} -mx $MAX_QUAL
		"""

rule assess_variant_plot:
	input:
		pr=list(pr_files)
	log:
		LOGS / "assess_mutref_calls/trimming_variant_summary.log"
	resources:
		mem="16GiB",
		runtime="30m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		figures=[
				FIGURES / f"assess/call/metrics/combo_variant_{metric}.png"
				for metric in ["f1", "recall", "precision"]
				],
		csv=TABLES / "assess/call/metrics/combo_variant_summary.csv"
	script:
		"../scripts/plot_variant_metrics.py"

rule assess_variant_average:
	input:
		csv=rules.assess_variant_plot.output.csv
	log:
		LOGS / "assess_mutref_calls/trimming_variant_summary_averages.log"
	resources:
		mem="2GiB",
		runtime="5m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		csv=TABLES / "assess/call/metrics/combo_variant_summary_averages.csv"
	script:
		"../scripts/average_variant_metrics.py"

rule assess_variant_fnfp:
	input:
		pr=list(pr_files)
	log:
		LOGS / "assess_mutref_calls/trimming_variant_fnfp.log"
	resources:
		mem="16GiB",
		runtime="15m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		csv=TABLES / "assess/call/metrics/combo_variant_fnfp.csv",
		fn_plot=FIGURES / "assess/call/metrics/combo_variant_delta_fn.png",
		fp_plot=FIGURES / "assess/call/metrics/combo_variant_delta_fp.png"
	script:
		"../scripts/extract_plot_fnfp.py"

rule assess_assembly_quast:
	input:
		assembly=rules.assembly_flye.output.assembly,
		reference=get_reference_genome
	log:
		LOGS / "assess/assembly/{tool}-{trimmer}/quast/{depth}x/{model}/{sample}.{tool}-{trimmer}.log"
	threads: 4
	resources:
		mem="32GiB",
		runtime="10m"
	params:
		one2one="--ambiguity-usage=one" # explicitly set this, although this option is default
	conda:
		ENVS / "quast.yaml"
	output:
		report_tsv= RESULTS / "assess/assembly/quast/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.report.tsv",  
		report_html= RESULTS / "assess/assembly/quast/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.report.html",
		icarus_html= RESULTS / "assess/assembly/quast/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.icarus.html",  
		icarus_helper1= RESULTS / "assess/assembly/quast/{tool}-{trimmer}/{depth}x/{model}/{sample}/icarus_viewers/contig_size_viewer.html", 
		icarus_helper2= RESULTS / "assess/assembly/quast/{tool}-{trimmer}/{depth}x/{model}/{sample}/icarus_viewers/alignment_viewer.html"
	shell:
		"""
		tmp_dir=$(mktemp -d)

		quast --threads {threads} {params.one2one} -r {input.reference} --output-dir $tmp_dir {input.assembly} 2> {log}
		mv ${{tmp_dir}}/report.tsv {output.report_tsv} 2>> {log}
		mv ${{tmp_dir}}/report.html {output.report_html} 2>> {log} 
		mv ${{tmp_dir}}/icarus.html {output.icarus_html} 2>> {log} 
		mv ${{tmp_dir}}/icarus_viewers/contig_size_viewer.html {output.icarus_helper1} 2>> {log} 
		mv ${{tmp_dir}}/icarus_viewers/alignment_viewer.html {output.icarus_helper2} 2>> {log}

		rm -rf $tmp_dir
		"""

# I had an issue installing sassy-rs using conda directive in snakemake rule,
# so i installed manually in my snakemake conda env
# pip install biopython sassy-rs --no-cache-dir --force-reinstall
rule assess_assembly_contam:
	input:
		assembly=rules.assembly_flye.output.assembly,
		info=rules.assembly_flye.output.info,
		contaminants=DATA / "contaminants"
	log:
		LOGS / "assess/assembly/contaminant/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.contaminants.log"
	threads: 4
	resources:
		mem="64GiB",
		runtime="1h"
	params:
		min_cov=0.90,
		min_id=0.90,
		kit=lambda wildcards: get_sequencing_kits(wildcards) 
	conda:
		ENVS / "align.yaml"
	output:
		contaminants=RESULTS / "assess/assembly/contaminant/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.contaminants.tsv"
	script:
		"../scripts/detect_adapter_contamination.py"

rule aggregate_assembly_contam:
	input:
		contam_report=list(contam_report_files)
	log:
		LOGS / "assess/assembly/contaminant/aggregate_contaminants.log"
	threads: 4
	resources:
		mem="8GiB",
		runtime="20m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		summary=TABLES / "assess/assembly/metrics/combo_contaminant_summary_count.csv",
		details=TABLES / "assess/assembly/metrics/combo_contaminant_details.csv"
	script:
		"../scripts/aggregate_contaminants.py"

rule plot_assembly_quast_metrics:
	input:
		reports=expand(
				RESULTS / "assess/assembly/quast/{combo}/{depth}x/{model}/{sample}.{combo}.report.tsv",
				combo=COMBINATIONS,
				depth=DEPTHS,
				sample=SAMPLES,
				model=MODELS)
	log:
		LOGS / "assess/assembly/plot_quast_metrics.log"
	resources:
		mem="16GiB",
		runtime="20m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		csv=TABLES / "assess/assembly/metrics/combo_quast_compiled_metrics.csv",
		error_plot=FIGURES / "assess/assembly/metrics/combo_assembly_errors_per_100kbp.png",
		nga50_plot=FIGURES / "assess/assembly/metrics/combo_assembly_nga50.png",
		delta_error_plot=FIGURES / "assess/assembly/metrics/combo_assembly_delta_errors.png",
		delta_nga50_plot=FIGURES / "assess/assembly/metrics/combo_assembly_delta_nga50.png"
	script:
		"../scripts/plot_assembly_metrics.py"
