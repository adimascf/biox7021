REPEAT = int(config.get("repeat", 1))

rule assembly_flye:
	input:
		reads=rules.downsample_rasusa.output.reads
	log:
		LOGS / "assembly/{trimmer}/{depth}x/{sample}.{trimmer}.flye.log"
	threads: 4
	resources:
		mem="128GiB",
		runtime=f"{5 * REPEAT}h"
	conda:
		ENVS / "flye.yaml"
	params:
		ont="--nano-hq",
	output:
		assembly=RESULTS / "assembly/{trimmer}/{depth}x/{sample}/{sample}.{trimmer}.{depth}x.assembly.fasta",
		graph=RESULTS / "assembly/{trimmer}/{depth}x/{sample}/{sample}.{trimmer}.{depth}x.assembly_graph.gfa",
		info=RESULTS / "assembly/{trimmer}/{depth}x/{sample}/{sample}.{trimmer}.{depth}x.assembly_info.txt"
	benchmark:
		repeat(BENCHMARK / "assembly/{trimmer}/{depth}x/{sample}.{trimmer}.tsv", REPEAT)
	shell:
		"""
		tmp_results=$(mktemp -d)
		flye {params.ont} {input.reads} --out-dir $tmp_results --threads {threads} 2> {log}
		mv "${{tmp_results}}"/assembly.fasta {output.assembly} 2>> {log}
		mv "${{tmp_results}}"/assembly_graph.gfa {output.graph} 2>> {log}
		mv "${{tmp_results}}"/assembly_info.txt {output.info} 2>> {log}

		rm -rf $tmp_results
		"""
	
