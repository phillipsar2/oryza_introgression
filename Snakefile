import config

### Wildcards ###
SAMPLE = glob_wildcards("data/samples/final_set/{sample}_1.fq.gz").sample
#print(SAMPLE)

REF = ["GLUM"]

# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

rule all:
    input:
#        fastqc = directory("reports/fastqc")
#        html=expand("qc/fastqc/{fastq}.html", fastq = FASTQ),
#        zip=expand("qc/fastqc/{fastq}_fastqc.zip", fastq = FASTQ),
        # Trim reads
#        trim = expand("data/rawdata/trimmed/{sample}.trim_1.fq.gz", sample = SAMPLE),
        # Aligning reads
        markdups = expand(config.mark_dups_outbam, sample = SAMPLE, REF=REF),
        # Assess quality of mapped reads
        bamqc = expand("reports/bamqc/{sample}_{REF}_stats/qualimapReport.html", sample = SAMPLE, REF=REF),
        # SNP calling
#        hap_caller = expand("data/gvcf/{REF}/{sample}.{REF}.g.vcf", sample = SAMPLE, REF=REF),
#        split_int = expand("data/processed/scattered_intervals/{REF}/{count}-scattered.interval_list", count = INTERVALS, REF=REF)
#        combine = expand(directory("data/interim/combined_database_bpres/{REF}/{count}"), count = INTERVALS, REF=REF),
#        joint_geno = expand("data/raw/vcf_bpres/{REF}/{count}.raw.{REF}.vcf", count = INTERVALS, REF=REF),
        # Process VCFs
#        get_snps = expand("data/raw/vcf_bpres/{REF}/{count}.raw.snps.{REF}.vcf", count = INTERVALS, REF=REF),
#        hard_filter = expand("data/processed/filtered_snps_bpres/{REF}/{count}.filtered.snps.{REF}.vcf", count = INTERVALS, REF=REF),
#        diagnostics = expand("reports/filtering/gvcf_{count}.{REF}.table", count = INTERVALS, REF=REF)
#        dp = expand("data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp1.snps.{REF}.vcf", count = INTERVALS, REF=REF),
#        dp2 = expand("data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.snps.{REF}.vcf", count = INTERVALS, REF=REF),
#        bgzip_vcf = expand("data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.snps.{REF}.vcf.gz", count = INTERVALS, REF=REF),
#        combine = expand("data/processed/filtered_snps_bpres/{REF}/oryza.{REF}.vcf.gz", REF=REF)
###        depth_diag = "reports/filtering_bpres/oryza_glum.table"
        # Combining fitlering SNPs and INDELs all together - not just the SNPs - then merging the files together
#        wholegenome = expand("data/processed/filtered_snps_bpres/{REF}/wholegenome.{REF}.vcf.gz", REF=REF)

# Rules
include: "rules/mapping.smk"
include: "rules/process_bam.smk"
#include: "rules/calling.smk"
#include: "rules/filtering.smk"
