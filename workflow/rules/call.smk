rule align_self:
	input:
		reads=rules.downsample_rasusa.output.reads,
		reference=rules.faidx_self.output.fasta,
	log:
		LOGS / "alignment/self/{trimmer}/{depth}x/{sample}.{trimmer}.minimap2.log"
	threads: 4
	resources:
		mem="32GiB",
		runtime="5m"
	conda:
		ENVS / "align.yaml"
	params:
		preset="map-ont",
		opts="-aL --cs --MD"
	output:
		bam=temp(RESULTS / "alignment/self/{trimmer}/{depth}x/{sample}.{trimmer}.minimap2.sorted.bam"),
		bai=temp(RESULTS / "alignment/self/{trimmer}/{depth}x/{sample}.{trimmer}.minimap2.sorted.bam.bai")
	shell:
		"""
		(minimap2 {params.opts} -t {threads} -x {params.preset} {input.reference} {input.reads} | \
				samtools sort -@ {threads} -o {output.bam}) 2> {log}
		samtools index {output.bam} -o {output.bai} 2>> {log}
		"""

use rule align_self as align_mutref with:
	input:
		reads=rules.downsample_rasusa.output.reads,
		reference=rules.faidx_mutref.output.fasta
	log:
		LOGS / "alignment/mutref/{trimmer}/{depth}x/{sample}.{trimmer}.minimap2.log"
	output:
		bam=temp(RESULTS / "alignment/mutref/{trimmer}/{depth}x/{sample}.{trimmer}.minimap2.sorted.bam"),
		bai=temp(RESULTS / "alignment/mutref/{trimmer}/{depth}x/{sample}.{trimmer}.minimap2.sorted.bam.bai")

rule call_self:
	input:
		alignment=rules.align_self.output.bam,
		bai=rules.align_self.output.bai,
		reference=rules.align_self.input.reference,
		faidx=rules.faidx_self.output.faidx
	log:
		LOGS / "calling/self/{trimmer}/{depth}x/{sample}.{trimmer}.clair3.log"
	threads: 4
	resources:
		mem="128GiB",
		runtime="6h"
	container:
		"docker://hkubal/clair3:v2.0.0"
	shadow:
		"shallow"
	output:
		vcf=RESULTS / "calling/self/{trimmer}/{depth}x/{sample}.{trimmer}.clair3.vcf.gz"
	script:
		"../scripts/calling/clair3.sh"

use rule call_self as call_mutref with:
	input:
		alignment=rules.align_mutref.output.bam,
		bai=rules.align_mutref.output.bai,
		reference=rules.align_mutref.input.reference,
		faidx=rules.faidx_mutref.output.faidx
	log:
		LOGS / "calling/mutref/{trimmer}/{depth}x/{sample}.{trimmer}.clair3.log"
	output:
		vcf=RESULTS / "calling/mutref/{trimmer}/{depth}x/{sample}.{trimmer}.clair3.vcf.gz"
