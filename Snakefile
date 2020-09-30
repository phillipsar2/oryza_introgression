import config

# fastq sample names for fastqc
#FASTQ = glob_wildcards("data/rawdata/{fastq}_fq.gz").fastq

# Sample names
#SAMPLE = glob_wildcards("data/rawdata/{sample}-R1.fastq").sample
SAMPLE = glob_wildcards("data/rawdata/trimmed/{sample}.trim_1.fq.gz").sample
# Testing a singular sample
#SAMPLE = ["IRIS_313-10274_120220_I297_FCC0DJBACXX_L6_RICwdsRSYHSD24-3-IPAAPEK-84"]
#print(SAMPLE)

# Maf names
MAF = glob_wildcards("data/ensemble_alignments/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.maf").maf
#print(MAF)

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
#        markdups = expand("data/interm/mark_dups/trim/{sample}.dedup.bam", sample = SAMPLE),
        # Assess quality of mapped reads
#        bamqc = expand("reports/bamqc/fastp/{sample}_stats/qualimapReport.html", sample = SAMPLE),
        # SNP calling
#        hap_caller = expand(config.haplo_caller, sample = SAMPLE),
#        split_int = expand(config.split_intervals, count = INTERVALS)
#        combine = expand(directory("data/interim/combined_database_bpres/domesticated/trim/{count}"), count = INTERVALS),
#        joint_geno = expand(config.joint_out, count = INTERVALS),
#        get_snps = expand(config.get_snps_out, count = INTERVALS),
#        hard_filter = expand(config.hard_filt, count = INTERVALS),
#        diagnostics = expand(config.diag, count = INTERVALS)
#        dp = expand(config.filt_dp1, count = INTERVALS),
#        dp2 = expand(config.filt_dp2, count = INTERVALS),
#        bgzip_vcf = expand(config.bgz_vcf, count = INTERVALS),
#        combine = config.combine,
#        depth_diag = "reports/filtering_bpres/oryza_glum.table"
        # Combining fitlering SNPs and INDELs all together - not just the SNPs - then merging the files together
#        wholegenome = config.whole_genome,
        # Liftover 
#        swap_maf = expand("data/ensemble_alignments/swap/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.swap.maf", maf = MAF),
#        convert__maf = expand("data/ensemble_alignments/chain/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.chain", maf = MAF),
#        combine_chain = "data/ensemble_alignments/merged_chains/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.combined.chain",
#        rename_vcfs = "data/vcf/oryza_glum.12chrom.vcf.gz",
        liftover_vcf = "data/lifted_vcf/oryza_sativa.half.lifted.vcf",

# Rules
#include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/calling.smk"
#include: "rules/filtering.smk"
include: "rules/liftover.smk"
