rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = "data/interm/mark_dups/{sample}.dedup.bam"
    output:
        outdir = "data/gvcf/{sample}.g.vcf"
    params:
        regions = "data/genome/contig.list"
    run:
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output.outdir} \
        --reference {input.ref} \
        --G StandardAnnotation \
        -G AS_StandardAnnotation \
        -L {params.regions} \
        -ERC BP_RESOLUTION")


# Scatter reference into intervals using SplitIntervals. Set intervals in Snakefile.
# https://gatk.broadinstitute.org/hc/en-us/articles/360036348372-SplitIntervals
# ``--scatter-count n` splits reference into n intervals
## This rule doesn't work right now
#rule split_intervals:
#    input:
#        ref = config.ref
#    output:
#        expand("data/processed/scattered_intervals/{count}-scattered.intervals", count = INTERVALS)
#    params:
#        regions = "data/genome/contig.list",
#        count = len(INTERVALS),
#        outdir = "data/processed/scattered_intervals"
#    run:
#        shell("gatk SplitIntervals \
#        -R {input.ref} \
#        -L {params.regions} \
#        --scatter-count {params.count} \
#        -O {params.outdir}")

# Combine GVCFs with GenomicsDBImport
# sample_map file is samplename with path/to/gvcf in the following line
# scattered interval file created by SpitIntevals

rule combine_gvcfs:
    input:
        gvcfs = expand("data/gvcf/{sample}.g.vcf", sample = SAMPLE),
        region = "data/processed/scattered_intervals/{count}-scattered.interval_list",
        map = "data/processed/sample_map"
    output:
        directory("data/interim/combined_database_bpres/{count}")
    params:
#        region = "data/processed/scattered_intervals/{count}-scattered.intervals_list",
        tmp = "/scratch/aphillip/genomicsdbimport/{count}"
#        map = "data/processed/sample_map"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk --java-options \"-Xmx80g -Xms80g\" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output} \
        --batch-size 50 \
        --reader-threads 8 \
        --sample-name-map {input.map} \
        --intervals {input.region} \
        --tmp-dir {params.tmp}")
        shell("rm -rf {params.tmp}")


rule joint_geno:
    input:
        dir = directory("data/interim/combined_database_bpres/{count}"),
        ref = config.ref
    output:
        "data/raw/vcf_bpres/{count}.raw.vcf"
    params:
        db = "gendb://data/interim/combined_database_bpres/{count}",
        region = "data/processed/scattered_intervals/{count}-scattered.interval_list"
    run:
        shell("gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {params.db} \
        -L {params.region} \
        -new-qual \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --include-non-variant-sites \
        -O {output}")
