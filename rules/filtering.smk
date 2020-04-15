# Extract SNPs

rule get_snps:
    input:
        ref = config.ref,
        vcf = "data/raw/vcf_bpres/{count}.raw.vcf"
    output:
        "data/raw/vcf_bpres/{count}.raw.snps.vcf"
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
#        vcf = ancient("data/raw/vcf_bpres/{count}.raw.snps.vcf")        
        vcf = "data/raw/vcf_bpres/{count}.raw.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{count}.filtered.snps.vcf"
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
        vcf = "data/processed/filtered_snps_bpres/{count}.filtered.snps.vcf",
        ref = config.ref
    output:
        "reports/filtering/gvcf_{count}.table"
    run:
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -O {output}")


# Filter SNPs to only biallelic and sites shared by all induviduals

rule filter_nocall:
    input:
        ref = config.ref,
        vcf = "data/processed/filtered_snps_bpres/{count}.filtered.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{count}.filtered.nocall.vcf"
    run: 
        shell("gatk SelectVariants -V {input.vcf} --max-nocall-fraction 0 --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}")
        

# Diagnositc filters look good. Continue filtering for depth. 

rule filter_depth:
    input:
        vcf = "data/processed/filtered_snps_bpres/{count}.filtered.nocall.vcf",
        ref = config.ref
    output:
        dp = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.snps.vcf",
        dp2 = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < 3 || DP > 77\" \
        -G-filter-name \"DP_3-77\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}")
        shell("gatk SelectVariants -V {output.dp} --max-nocall-fraction 0 --exclude-filtered true --restrict-alleles-to BIALLELIC -O {output.dp2}")



# Use vcftools to filter max DP less than or equal to mean DP (for each induvidual)
# Not in use

rule filter_maxdepth:
    input:
        vcf = "data/processed/filtered_snps_bpres/{count}.filtered.dp_min3.snps.vcf",
#        dp = "data/processed/filtered_snps_bpres/{count}.filtered.dp_min3_maxmean.snps.vcf"
    output:
        dp = "data/processed/filtered_snps_bpres/{count}.filtered.dp_min3_maxmean.snps.vcf"
#        dp2 = "data/processed/filtered_snps_bpres/{count}.filtered.dp_min3_maxmean.nocall.snps.vcf"  
    params:
        out = "{count}.filtered.dp_min3_maxmean.snps"
    run:
        shell("vcftools --vcf {input} --max-meanDP -c > {output.dp}")
#        shell("gatk SelectVariants -V {input.dp} --max-nocall-fraction 0 --exclude-filtered true --restrict-alleles-to BIALLELIC -O {output.dp2}")


rule bgzip_vcf:
    input:
        "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf.gz"
    run:
        shell("bgzip {input}")
        shell("tabix -p vcf {output}")


rule combine_vcfs:
    input:
        expand("data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf.gz", count = INTERVALS)
    output:
        "data/processed/filtered_snps_bpres/oryza_glum.vcf.gz"
    run:
        shell("bcftools concat {input} -Oz -o {output}")


rule depth:
    input:
#        vcf = "data/processed/filtered_snps_bpres/oryza_glum.vcf.gz",
        vcf = "data/processed/filtered_snps_bpres/depth_analysis.vcf.gz",
        ref = config.ref
    output:
#        dp = "reports/filtering_bpres/oryza_glum.table"
        dp = "reports/filtering_bpres/depth_analysis.table"
    run:
        shell("tabix -p vcf {input.vcf}")
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -GF DP \
        -O {output.dp}")


# Filter the whole-genome file

rule filter_wholegenome:
    input:
        vcf = "data/raw/vcf_bpres/{count}.raw.vcf",
        ref = config.ref
    output:
        dp = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.wholegenome.vcf",
        dp2 = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.wholegenome.vcf"
    run:
        shell("gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < 3 || DP > 77\" \
        -G-filter-name \"DP_3-77\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}")
        shell("gatk SelectVariants -V {output.dp} --max-nocall-fraction 0 --exclude-filtered true --restrict-alleles-to BIALLELIC -O {output.dp2}")

rule combine_wgenomevcfs:
    input:
        expand("data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.wholegenome.vcf", count = INTERVALS)
    output:
        "data/processed/filtered_snps_bpres/oglum_wholegenome.vcf.gz"
    run:
        shell("bcftools concat {input} -Oz -o {output}")
