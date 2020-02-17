import os

configfile: "gombe_config.json"

fastq_directory = "/mnt/storage/SAYRES/GOMBE/fastqs"
temp_directory = "temp"

fastqc_path = "fastqc"
multiqc_path = "multiqc"
bbduksh_path = "bbduk.sh"
samtools_path = "samtools"
bwa_path = "bwa"
gatk_path = "gatk"
picard_path = "picard"
sambamba_path = "sambamba"
bedtools_path = "bedtools"
sort_bed_path = "sort-bed"
bcftools_path = "bcftools"
tabix_path = "tabix"

fastq_prefixes = [
	config[x]["fq1"][:-9] for x in config["file_ids"]] + [
		config[x]["fq2"][:-9] for x in config["file_ids"]]

pairwise_list = [
	"s7507_calculus_exome-s7507_dentine_exome",
	"s7365_feces_exome-s7365_urine_exome",
	"s7150_feces_exome-s7150_urine_exome",
	"s7507_feces_exome-s7507_urine_exome"]

rule all:
	input:
		"multiqc/multiqc_report.html",
		"trimmed_multiqc/multiqc_report.html",
		expand(
			"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai",
			sample=config["final_id_list"], genome=["pantro4"]),
		expand(
			"stats/{sample}.{genome}.sorted.merged.mkdup.bam.stats",
			sample=config["final_id_list"], genome=["pantro4"]),
		expand(
			"stats/{sample}.{genome}.sorted.merged.mkdup.bam.NODUPS.stats",
			sample=config["final_id_list"], genome=["pantro4"]),
		expand(
			"results/{sample}.{genome}.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
			sample=config["final_id_list"], genome=["pantro4"]),
		expand(
			"vcf_compare/results.{comparison}.{genome}.{depth}.vcfcompare.norandom_or_unplaced.csv",
			comparison=pairwise_list, genome=["pantro4"], depth=["4", "6", "8", "10"]),
		# add downsampled
		expand(
			"downsampled_stats/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.stats",
			sample=config["to_downsample"], genome=["pantro4"]),
		expand(
			"downsampled_stats/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.NODUPS.stats",
			sample=config["to_downsample"], genome=["pantro4"]),
		expand(
			"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
			sample=config["to_downsample"], genome=["pantro4"])

rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		shell(
			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell(
			"{params.samtools} faidx {output.new}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {output.new}")
		# bwa
		shell(
			"{params.bwa} index {output.new}")

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{prefix}.fastq.gz")
	output:
		"fastqc/{prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc fastqc"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe qtrim=rl trimq=10 minlength=30"

rule fastqc_analysis_trimmed_paired:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "trimmed_fastqc/{sample}_trimmed_read1_fastqc.html",
		html2 = "trimmed_fastqc/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o trimmed_fastqc {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"trimmed_fastqc/{sample}_trimmed_{read}_fastqc.html",
			sample=config["file_ids"], read=["read1", "read2"])
	output:
		"trimmed_multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o trimmed_multiqc trimmed_fastqc"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{genome}.fasta",
		fai = "reference/{genome}.fasta.fai",
	output:
		"processed_bams/{sample}.{genome}.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	shell:
		" {params.bwa} mem -t {params.threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule merge_bams:
	input:
		bams = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.sorted.bam",
			sample=config["to_merge"][wildcards.sample], genome=wildcards.genome),
		bais = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.sorted.bam.bai",
			sample=config["to_merge"][wildcards.sample], genome=wildcards.genome)
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	threads:
		4
	params:
		sambamba = sambamba_path,
		threads = 4
	run:
		if len(input.bams) > 1:
			shell("{params.sambamba} merge -t {params.threads} {output} {input.bams}")
		else:
			shell(
				"ln -s ../{} {{output}} && touch -h {{output}}".format(input.bams[0]))

rule index_merged_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		metrics = "results/{sample}.{genome}.picard_mkdup_metrics.txt"
	threads: 4
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.merged.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule bam_stats_nodups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.merged.mkdup.bam.NODUPS.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} view -F 1024 -b {input.bam} | "
		"{params.samtools} stats - | grep ^SN | cut -f 2- > {output}"

rule gatk_haplotype_caller:
	input:
		ref = "reference/{genome}.fasta",
		fai = "reference/{genome}.fasta.fai",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"vcfs/{sample}.{genome}.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		8
	shell:
		"""{params.gatk} --java-options "-Xmx30g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -O {output}"""

rule bcftools_merge_vcfs:
	input:
		vcf1 = "vcfs/{comp1}.{genome}.raw.vcf.gz",
		vcf2 = "vcfs/{comp2}.{genome}.raw.vcf.gz"
	output:
		"merged_vcfs/{comp1}-{comp2}.{genome}.raw.vcf.gz"
	params:
		bcftools = bcftools_path
	shell:
		"{params.bcftools} merge -m all -o {output} -O z {input}"

rule index_merged_vcf:
	input:
		"merged_vcfs/{comp1}-{comp2}.{genome}.raw.vcf.gz"
	output:
		"merged_vcfs/{comp1}-{comp2}.{genome}.raw.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule gatk_haplotype_caller_second_pass:
	input:
		ref = "reference/{genome}.fasta",
		fai = "reference/{genome}.fasta.fai",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai",
		intervals = "merged_vcfs/{comp1}-{comp2}.{genome}.raw.vcf.gz",
		interval_tbi = "merged_vcfs/{comp1}-{comp2}.{genome}.raw.vcf.gz.tbi"
	output:
		"final_vcfs/{sample}.{genome}.intervals-{comp1}-{comp2}.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		8
	shell:
		"""{params.gatk} --java-options "-Xmx30g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} """
		"""-ERC BP_RESOLUTION -L {input.intervals}"""

rule filter_genome_cov_depth4:
	input:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov"
	output:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov.dp4.bed"
	shell:
		"""awk '($4 >= 4)' {input} > {output}"""

rule filter_genome_cov_depth8:
	input:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov"
	output:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov.dp8.bed"
	shell:
		"""awk '($4 >= 8)' {input} > {output}"""

rule intersect_genome_cov:
	input:
		bed_a = "results/{comp1}.{genome}.mapq20_noDup.genome_cov.{depth}.bed",
		bed_b = "results/{comp2}.{genome}.mapq20_noDup.genome_cov.{depth}.bed"
	output:
		"results/{comp1}-{comp2}.{genome}.{depth}.bed"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} intersect -a {input.bed_a} -b {input.bed_b} > {output}"

# rule compare_vcfs_using_refpanel:
# 	input:
# 		vcf_a = "final_vcfs/{comp1}.{genome}.intervals-{comp1}-{comp2}.raw.vcf.gz",
# 		vcf_b = "final_vcfs/{comp2}.{genome}.intervals-{comp1}-{comp2}.raw.vcf.gz",
# 		vcf_b_idx = "vcfs/{comp2}.{genome}.filtered.sorted.vcf.gz.tbi",
# 		bed_a = "results/{comp1}.{genome}.mapq20_noDup.genome_cov.{depth}.bed",
# 		bed_b = "results/{comp2}.{genome}.mapq20_noDup.genome_cov.{depth}.bed",
# 		ref_panel = ref_panel_vcf
# 	output:
# 		"vcf_compare/{comp1}-{comp2}.{genome}.{depth}.vcfcompare.table.csv"
# 	shell:
# 		"python scripts/Compare_vcfs_with_refpanel_and_beds.py "
# 		"--vcf1 {input.vcf_a} "
# 		"--vcf2 {input.vcf_b} "
# 		"--bed1 {input.bed_a} "
# 		"--bed2 {input.bed_b} "
# 		"--vcf_ref_panel {input.ref_panel} "
# 		"--output_file {output}"

rule compare_vcfs:
	input:
		vcf_a = "final_vcfs/{comp1}.{genome}.intervals-{comp1}-{comp2}.raw.vcf.gz",
		vcf_b = "final_vcfs/{comp2}.{genome}.intervals-{comp1}-{comp2}.raw.vcf.gz",
		vcf_a_idx = "vcfs/{comp1}.{genome}.filtered.sorted.vcf.gz.tbi",
		vcf_b_idx = "vcfs/{comp2}.{genome}.filtered.sorted.vcf.gz.tbi"
	output:
		"vcf_compare/{comp1}-{comp2}.{genome}.vcfcompare.table.csv"
	shell:
		"python scripts/Compare_vcfs.py "
		"--vcf1 {input.vcf_a} "
		"--vcf2 {input.vcf_b} "
		"--output_file {output}"

rule remove_unplaced_scaffolds:
	input:
		"vcf_compare/{comp1}-{comp2}.{genome}.vcfcompare.table.csv"
	output:
		"vcf_compare/{comp1}-{comp2}.{genome}.vcfcompare.table.norandom_or_unplaced.csv"
	shell:
		"grep -v '_random' {input} | grep -v 'chrUn' > {output}"

rule process_dropout:
	input:
		"vcf_compare/{comp1}-{comp2}.{genome}.vcfcompare.table.norandom_or_unplaced.csv"
	output:
		"vcf_compare/results.{comp1}-{comp2}.{genome}.{depth}.vcfcompare.norandom_or_unplaced.csv"
	params:
		one = "{comp1}",
		two = "{comp2}",
		dp = "{depth}",
		mq = "30",
		gq = "30"
	shell:
		"python scripts/Process_dropout.py --input {input} "
		"--dp {params.dp} "
		"--mq {params.mq} "
		"--gq {params.gq} "
		"--sample1_name {params.one} "
		"--sample2_name {params.two} "
		"--output {output}"

rule make_bedtools_genome_file:
	input:
		fai = "reference/{genome}.fasta.fai"
	output:
		"reference/{genome}.genome"
	shell:
		"awk -v OFS='\t' {{'print $1,$2'}} {input.fai} > {output}"

rule genome_cov:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		idx = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai",
		genome = "reference/{genome}.genome"
	output:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov"
	params:
		samtools = samtools_path,
		bedtools = bedtools_path
	shell:
		"{params.samtools} view {input.bam} -b -F 1024 -q 20 | "
		"{params.bedtools} genomecov -bg -ibam stdin -g {input.genome} > {output}"

rule bedops_sort_genome_cov:
	input:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov"
	output:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov.bedopssorted.bed"
	params:
		sort_bed = sort_bed_path
	shell:
		"{params.sort_bed} {input} > {output}"

rule bedtools_intersect_cds:
	input:
		sample = "results/{sample}.{genome}.mapq20_noDup.genome_cov.bedopssorted.bed",
		region = "misc/{genome}.ensembl.cds.bed"
	output:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov.bedopssorted.cds.INTERSECTION.bed"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} intersect -a {input.sample} -b {input.region} > {output}"

rule compute_histogram_from_bed:
	input:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov.bedopssorted.cds.INTERSECTION.bed"
	output:
		"results/{sample}.{genome}.mapq20_noDup.genome_cov.bedopssorted.cds.hist"
	shell:
		"python scripts/Compute_histogram_from_bed.py --bed {input} --outfile {output}"

rule downsample_bams:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		"downsampled/{sample}.{genome}.sorted.merged.downsampled.bam"
	params:
		samtools = samtools_path,
		downsample_fraction = lambda wildcards: config["downsample_fraction"][wildcards.sample]
	shell:
		"{params.samtools} view -s 0.{params.downsample_fraction} -b {input.bam} > {output}"

rule index_downsampled_merged_bam:
	input:
		"downsampled/{sample}.{genome}.sorted.merged.downsampled.bam"
	output:
		"downsampled/{sample}.{genome}.sorted.merged.downsampled.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups_downsampled:
	input:
		bam = "downsampled/{sample}.{genome}.sorted.merged.downsampled.bam",
		bai = "downsampled/{sample}.{genome}.sorted.merged.downsampled.bam.bai"
	output:
		bam = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam",
		metrics = "results/{sample}.{genome}.picard_mkdup_metrics.downsampled.txt"
	threads: 4
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam_downsampled:
	input:
		"downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam"
	output:
		"downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats_downsampled:
	input:
		bam = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam",
		bai = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.bai"
	output:
		"downsampled_stats/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule bam_stats_nodups_downsampled:
	input:
		bam = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam",
		bai = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.bai"
	output:
		"downsampled_stats/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.NODUPS.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} view -F 1024 -b {input.bam} | "
		"{params.samtools} stats - | grep ^SN | cut -f 2- > {output}"

rule genome_cov_downsampled:
	input:
		bam = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam",
		idx = "downsampled/{sample}.{genome}.sorted.merged.downsampled.mkdup.bam.bai",
		genome = "reference/{genome}.genome"
	output:
		"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov"
	params:
		samtools = samtools_path,
		bedtools = bedtools_path
	shell:
		"{params.samtools} view {input.bam} -b -F 1024 -q 20 | "
		"{params.bedtools} genomecov -bg -ibam stdin -g {input.genome} > {output}"

rule bedops_sort_genome_cov_downsampled:
	input:
		"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov"
	output:
		"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov.bedopssorted.bed"
	params:
		sort_bed = sort_bed_path
	shell:
		"{params.sort_bed} {input} > {output}"

rule bedtools_intersect_cds_downsampled:
	input:
		sample = "downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov.bedopssorted.bed",
		region = "misc/{genome}.ensembl.cds.bed"
	output:
		"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.INTERSECTION.bed"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} intersect -a {input.sample} -b {input.region} > {output}"

rule compute_histogram_from_bed_downsampled:
	input:
		"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.INTERSECTION.bed"
	output:
		"downsampled_results/{sample}.{genome}.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist"
	shell:
		"python scripts/Compute_histogram_from_bed.py --bed {input} --outfile {output}"

# rule filter_individual_vcfs:
# 	input:
# 		"vcfs/{sample}.{genome}.raw.vcf.gz"
# 	output:
# 		"vcfs/{sample}.{genome}.filtered.vcf.gz"
# 	params:
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} filter -O z --include "
# 		"'INFO/MQ>=20 && QUAL>=30 && FORMAT/GQ>=30' {input} > {output}"
#
# rule index_filtered_vcf:
# 	input:
# 		"vcfs/{sample}.{genome}.filtered.vcf.gz"
# 	output:
# 		"vcfs/{sample}.{genome}.filtered.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"
#
# rule sort_filtered_vcf:
# 	input:
# 		vcf = "vcfs/{sample}.{genome}.filtered.vcf.gz",
# 		tbi = "vcfs/{sample}.{genome}.filtered.vcf.gz.tbi"
# 	output:
# 		"vcfs/{sample}.{genome}.filtered.sorted.vcf.gz"
# 	params:
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} sort -O z --output-file {output} {input.vcf}"
#
# rule index_sorted_vcf:
# 	input:
# 		"vcfs/{sample}.{genome}.filtered.sorted.vcf.gz"
# 	output:
# 		"vcfs/{sample}.{genome}.filtered.sorted.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"
