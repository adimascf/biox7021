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

rule assess_mutref_calls:
	input:
		query_vcf=RESULTS / "calling/mutref/{trimmer}/{sample}.{trimmer}.clair3.vcf.gz",
		truth_vcf=get_truth_vcf,
		mutreference=get_mutreference_genome,
		faidx=rules.faidx_mutref.output.faidx,
		bed=rules.make_full_bed.output.bed
	log:
		LOGS / "assess_mutref_calls/{trimmer}/{sample}.{trimmer}.log"
	resources:
		mem="32GiB",
		runtime="30m"
	params:
		opts="--largest-variant 50 --credit-threshold 1.0",
		prefix=lambda wildcards: RESULTS / f"assess/mutref/{wildcards.trimmer}/{wildcards.sample}/{wildcards.sample}"
	output:
		pr_summary=RESULTS / "assess/mutref/{trimmer}/{sample}/{sample}.precision-recall-summary.tsv",
		pr=RESULTS / "assess/mutref/{trimmer}/{sample}/{sample}.precision-recall.tsv",
		summary=RESULTS / "assess/mutref/{trimmer}/{sample}/{sample}.summary.vcf",
		query=RESULTS / "assess/mutref/{trimmer}/{sample}/{sample}.query.tsv",
		truth=RESULTS / "assess/mutref/{trimmer}/{sample}/{sample}.truth.tsv"
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

rule best_f1:
	input:
		pr=list(pr_files)
	log:
		LOGS / "visualisation/f1_trimming_only.log"
	resources:
		mem="16GiB",
		runtime="30m"
	conda:
		ENVS / "f1_score.yaml"
	output:
		figures=[
				FIGURES / f"assess/mutref/plots/trimming_{metric}.png"
				for metric in ["f1", "recall", "precision"]
				],
		csv=TABLES / "assess/mutref/plots/trimming_summary.csv"
	script:
		"../scripts/best_f1.py"


