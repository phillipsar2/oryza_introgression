rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = "data/interm/mark_dups/trim/{sample}.dedup.bam"
    output:
        outdir = "data/gvcf/trim/{sample}.g.vcf"
    params:
        regions = config.contig_list
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
rule split_intervals:
    input:
        ref = config.ref
    output:
        expand("data/processed/scattered_intervals/domesticated/{count}-scattered.interval_list", count = INTERVALS)
#        config.ref_int
#        expand("data/processed/scattered_intervals/{count}-scattered.interval_list", count = INTERVALS)
    params:
        regions = config.contig_list,
        dir = config.int_dir
    run:
        shell("gatk SplitIntervals -R {input.ref} -L {params.regions} --scatter-count 200 -O {params.dir}")

# Combine GVCFs with GenomicsDBImport
# sample_map file is samplename with path/to/gvcf in the following line
# scattered interval file created by SpitIntevals

rule combine_gvcfs:
    input:
        gvcfs = expand("data/gvcf/trim/{sample}.g.vcf", sample = SAMPLE),
        region = config.int_region,
#        region = "data/processed/scattered_intervals/{count}-scattered.interval_list",
        map = config.sample_map
    output:
#        directory("data/interim/combined_database_bpres/{count}")
        directory("data/interim/combined_database_bpres/domesticated/trim/{count}")
#        config.comb_dir
    params:
#        region = "data/processed/scattered_intervals/{count}-scattered.intervals_list",
        tmp = "/scratch/aphillip/genomicsdbimport/{count}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk --java-options \"-Xmx90g -Xms90g\" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output} \
        --batch-size 50 \
        --reader-threads 8 \
        --sample-name-map {input.map} \
        --intervals {input.region} --tmp-dir {params.tmp}")
        shell("rm -rf {params.tmp}")


rule joint_geno:
    input:
#        dir = directory("data/interim/combined_database_bpres/{count}"),
        dir = directory(config.comb_dir),
        ref = config.ref
    output:
        config.joint_out
#        "data/raw/vcf_bpres/domesticated/{count}.raw.vcf"
    params:
#        db = "gendb://data/interim/combined_database_bpres/domesticated/{count}",
        db = config.db,
#        region = "data/processed/scattered_intervals/{count}-scattered.interval_list"
        region = config.int_region
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
