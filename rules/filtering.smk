# Extract SNPs

rule get_snps:
    input:
        ref = config.ref,
        vcf = "data/raw/vcf_bpres/{REF}/{count}.raw.{REF}.vcf"
    output:
        "data/raw/vcf_bpres/{REF}/{count}.raw.snps.{REF}.vcf"
    run:
        shell("gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}")


# Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225


# Apply the base GATK filter

rule filter_snps:
    input:
        ref = config.ref,
        vcf = "data/raw/vcf_bpres/{REF}/{count}.raw.snps.{REF}.vcf"
    output:
        "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.snps.{REF}.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QD < 2.0\" --filter-name \"QD2\" \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"SOR > 3.0\" --filter-name \"SOR3\" \
        -filter \"FS > 60.0\" --filter-name \"FS60\" \
        -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
        -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
        -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
        -O {output}")

# Filtering diagnostics
# Extract variant quality scores
# https://evodify.com/gatk-in-non-model-organism/

rule diagnostics:
    input: 
        vcf = "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.snps.{REF}.vcf",
        ref = config.ref
    output:
        "reports/filtering/gvcf_{count}.{REF}.table"
    run:
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -O {output}")


# Filter SNPs to only biallelic
# --max-nocall-fraction is 0 for O. glum and 0.33 for O. sativa

rule filter_nocall:
    input:
        ref = config.ref,
        vcf = "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.snps.{REF}.vcf"
    output:
        "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.nocall.{REF}.vcf"
    run: 
        shell("gatk SelectVariants -V {input.vcf} --max-nocall-fraction 0.5 --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}")
        

# Diagnositc filters look good. Continue filtering for depth. 

# Evaluate depth across samples to set DP filter

#rule depth:
#    input:
#        vcf = "data/processed/filtered_snps_bpres/oryza_glum.vcf.gz",
#        vcf = "data/processed/filtered_snps_bpres/depth_analysis.vcf.gz",
#        ref = config.ref
#    output:
#        dp = config.depth_table
#    run:
#        shell("tabix -p vcf {input.vcf}")
#        shell("gatk VariantsToTable \
#        -R {input.ref} \
#        -V {input.vcf} \
#        -F CHROM -F POS -GF DP \
#        -O {output.dp}")

# Fitlter by depth of each individual
# O. glum samples filtered 3 < DP < 77 (inbreeding)
# O. sativa samples filtered 3 < DP < 90 (intbreeding)

rule filter_depth:
    input:
        vcf =  "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.nocall.{REF}.vcf",
        ref = config.ref
    output:
        dp = "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp1.snps.{REF}.vcf"
    run:
        shell("gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < 3 || DP > 90\" \
        -G-filter-name \"DP_3-90\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}")

rule depth_nocall:
    input:
        "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp1.snps.{REF}.vcf"
    output:
        "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.snps.{REF}.vcf"
    run:
        shell("gatk SelectVariants -V {input} --exclude-filtered true -O {output}")

# gzip vcfs

rule bgzip_vcf:
    input:
        "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.snps.{REF}.vcf"
    output:
        "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.snps.{REF}.vcf.gz"
    run:
        shell("bgzip {input}")
        shell("tabix -p vcf {output}")

# Combine individual vcfs with bcftools

rule combine_vcfs:
    input:
        expand("data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.snps.{REF}.vcf.gz", count = INTERVALS, REF=REF)
    output:
        "data/processed/filtered_snps_bpres/{REF}/oryza.{REF}.vcf.gz"
    run:
        shell("bcftools concat {input} -Oz -o {output}")



# Filter the whole-genome file
# Change the DP filter maximum depending on the genome being analyzed

rule filter_wholegenome:
    input:
        vcf = "data/raw/vcf_bpres/{REF}/{count}.raw.{REF}.vcf",
        ref = config.ref
    output:
        dp = "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp1.wholegenome.{REF}.vcf",
        dp2 = "data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.wholegenome.{REF}.vcf"
    run:
        shell("gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < 3 || DP > 90\" \
        -G-filter-name \"DP_3-90\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}")
        shell("gatk SelectVariants -V {output.dp} --exclude-filtered true -O {output.dp2}")

rule combine_wgenomevcfs:
    input:
        expand("data/processed/filtered_snps_bpres/{REF}/{count}.filtered.dp2.nocall.wholegenome.{REF}.vcf", count = INTERVALS, REF=REF)
    output:
        "data/processed/filtered_snps_bpres/{REF}/wholegenome.{REF}.vcf.gz"
    run:
        shell("bcftools concat {input} -Oz -o {output}")
