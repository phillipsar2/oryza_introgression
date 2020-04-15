import config

# fastq sample names for fastqc
#FASTQ = glob_wildcards("data/rawdata/{fastq}_fq.gz").fastq

# Sample names
#SAMPLE = glob_wildcards("data/rawdata/{sample}-R1.fastq").sample
SAMPLE = glob_wildcards("data/rawdata/{sample}_1.fq.gz").sample

# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

rule all:
    input:
#        fastqc = directory("reports/fastqc")
#        html=expand("qc/fastqc/{fastq}.html", fastq = FASTQ),
#        zip=expand("qc/fastqc/{fastq}_fastqc.zip", fastq = FASTQ),
        trim = expand("data/rawdata/trimmed/{sample}.forward.1.fq", sample = SAMPLE)
        # Aligning reads
#        markdups = expand("data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
        # Assess quality of mapped reads
#        bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLE),
        # SNP calling
#        hap_caller = expand("data/gvcf/{sample}.g.vcf", sample = SAMPLE),
#        split_int = expand("data/processed/scattered_intervals/{count}-scattered.intervals", count = INTERVALS)
#        directory("data/interim/combined_database_bpres/200")
#        expand("data/raw/vcf_bpres/{count}.raw.vcf", count = INTERVALS)
#        base_filter = expand("data/processed/filtered_snps_bpres/{count}.filtered.snps.vcf", count = INTERVALS),
#        diagnostics = expand("reports/filtering/gvcf_{count}.table", count = INTERVALS)
#        dp = expand("data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.snps.vcf", count = INTERVALS),
#        dp2 = expand("data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf", count = INTERVALS),
#        bgzip_vcf = expand("data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf.gz", count = INTERVALS),
#        combine = "data/processed/filtered_snps_bpres/oryza_glum.vcf.gz",
#        depth_diag = "reports/filtering_bpres/oryza_glum.table"
#        wholegenome = "data/processed/filtered_snps_bpres/oglum_wholegenome.vcf.gz"

# Rules
include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/calling.smk"
#include: "rules/filtering.smk"
