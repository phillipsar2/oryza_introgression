### Specify reference genome ###
import config
#REF = ["SAT"]
REF = ["GLUM"]


### Wildcards ###
SAMPLE = glob_wildcards("data/samples/final_set/{sample}_1.fq.gz").sample
#print(SAMPLE)


# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

# Identify list of chromosomes as wildcards with a little python magic
CHROM = [line.strip() for line in open(config.contig_list, 'r')] 

rule all:
    input:
#        fastqc = directory("reports/fastqc")
#        html=expand("qc/fastqc/{fastq}.html", fastq = FASTQ),
#        zip=expand("qc/fastqc/{fastq}_fastqc.zip", fastq = FASTQ),
        # Trim reads
#        trim = expand("data/rawdata/trimmed/{sample}.trim_1.fq.gz", sample = SAMPLE),
        # Aligning reads
#        markdups = expand("data/interm/mark_dups/{REF}/{sample}.dedup.{REF}.bam", sample = SAMPLE, REF=REF),
        # Assess quality of mapped reads
#        bamqc = expand("reports/bamqc/{sample}_{REF}_stats/qualimapReport.html", sample = SAMPLE, REF=REF),
        # SNP calling
#        hap_caller = expand("data/gvcf/{REF}/{sample}.{REF}.g.vcf", sample = SAMPLE, REF=REF),
#        split_int = expand("data/processed/scattered_intervals/{REF}/{interval}-scattered.interval_list", interval = INTERVALS, REF=REF),
#        combine = directory(expand("data/interm/combined_database_bpres/{REF}/{interval}", interval = INTERVALS, REF=REF)),
#        joint_geno = expand("data/raw/vcf_bpres/{REF}/{interval}.raw.{REF}.vcf", interval = INTERVALS, REF=REF),
        # Process VCFs
#        get_snps = expand("data/raw/vcf_bpres/{REF}/{interval}.raw.snps.{REF}.vcf", interval = INTERVALS, REF=REF),
#        hard_filter = expand("data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.snps.{REF}.vcf", interval = INTERVALS, REF=REF),
#        diagnostics = expand("reports/filtering/gvcf_{interval}.{REF}.table", interval = INTERVALS, REF=REF)
#        filt_notcall = expand("data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.nocall.{REF}.vcf", interval = INTERVALS, REF=REF),
        # Filter for Depth
#        depth_diag = expand("reports/filtering/depth.{interval}.filtered.nocall.{REF}.table", interval = INTERVALS, REF=REF), 
#        filt_dp = expand(config.filt_depth, interval = INTERVALS, REF=REF),
#        filt_nocall = expand(config.dp_nocall, interval = INTERVALS, REF=REF),
#        bgzip_vcf = expand(config.bgzip_vcf, interval = INTERVALS, REF=REF),
#        combine = expand("data/processed/filtered_snps_bpres/{REF}/oryza.{REF}.vcf.gz", REF=REF)
###        depth_diag = "reports/filtering_bpres/oryza_glum.table"
        # Combining fitlering SNPs and INDELs all together - not just the SNPs - then merging the files together
#        wholegenome = expand("data/processed/filtered_snps_bpres/{REF}/wholegenome.{REF}.vcf.gz", REF=REF)
        # Population strucutre
#         beagle = expand(config.chrom_beagle, REF=REF, chrom=CHROM),
#         all_beagle = expand(config.combined_beagle, REF=REF),
         pca = expand(config.pca, REF=REF)

# Rules
#include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/calling.smk"
#include: "rules/filtering.smk"
include: "rules/pop_structure.smk"
