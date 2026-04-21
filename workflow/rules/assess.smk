rule make_full_bed:
	input:
		faidx=rules.faidx_mutref.output.faidx
	log:
		LOGS / "make_fule_bed/{sample}.bed"
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
		query_vcf=RESULTS / "calling/{trimmer}/{depth}x/{sample}.{trimmer}.clair3.vcf.gz",
		truth_vcf=get_truth_vcf,  
		mutreference=get_mutreference_genome,
		faidx=rules.faidx_mutref.output.faidx,
		bed=rules.make_full_bed.output.bed
	log:
		LOGS / "assess_mutref_calls/{trimmer}/{depth}x/{sample}.{trimmer}.log"
	resources:
		mem="32GiB",
		runtime="30m"
	params:
		opts="--largest-variant 50 --credit-threshold 1.0",
		prefix=lambda wildcards: RESULTS / f"assess/call/{wildcards.trimmer}/{wildcards.depth}x/{wildcards.sample}/{wildcards.sample}"
	output:  
		pr_summary=RESULTS / "assess/call/{trimmer}/{depth}x/{sample}/{sample}.precision-recall-summary.tsv",
		pr=RESULTS / "assess/call/{trimmer}/{depth}x/{sample}/{sample}.precision-recall.tsv",
		summary=RESULTS / "assess/call/{trimmer}/{depth}x/{sample}/{sample}.summary.vcf",
		query=RESULTS / "assess/call/{trimmer}/{depth}x/{sample}/{sample}.query.tsv",
		truth=RESULTS / "assess/call/{trimmer}/{depth}x/{sample}/{sample}.truth.tsv"
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
		ENVS / "assess_variant_python.yaml"
	output:
		figures=[
				FIGURES / f"assess/call/metrics/trimming_variant_{metric}.png"
				for metric in ["f1", "recall", "precision"]
				],
		csv=TABLES / "assess/call/metrics/trimming_variant_summary.csv"
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
		ENVS / "assess_variant_python.yaml"
	output:
		csv=TABLES / "assess/call/metrics/trimming_variant_summary_averages.csv"
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
		ENVS / "assess_variant_python.yaml"
	output:
		csv=TABLES / "assess/call/metrics/trimming_variant_fnfp.csv"
	script:
		"../scripts/extract_fnfp_numbers.py"

rule assess_assembly_quast:
	input:
		assembly=rules.assembly_flye.output.assembly,
		reference=get_reference_genome
	log:
		LOGS / "assess/assembly/{trimmer}/{depth}x/{sample}.{trimmer}.log"
	threads: 4
	resources:
		mem="32GiB",
		runtime="10m"
	params:
		one2one="--ambiguity-usage=one" # explicitly set this, although this option is default
	conda:
		ENVS / "quast.yaml"
	output:
		report_tsv= RESULTS / "assess/assembly/{trimmer}/{depth}x/{sample}/{sample}.{trimmer}.report.tsv",  
		report_html= RESULTS / "assess/assembly/{trimmer}/{depth}x/{sample}/{sample}.{trimmer}.report.html",
		icarus_html= RESULTS / "assess/assembly/{trimmer}/{depth}x/{sample}/{sample}.{trimmer}.icarus.html",  
		icarus_helper1= RESULTS / "assess/assembly/{trimmer}/{depth}x/{sample}/icarus_viewers/contig_size_viewer.html", 
		icarus_helper2= RESULTS / "assess/assembly/{trimmer}/{depth}x/{sample}/icarus_viewers/alignment_viewer.html"
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

