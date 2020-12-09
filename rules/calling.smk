rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = "data/interm/mark_dups/{REF}/{sample}.dedup.{REF}.bam"
    output:
        outdir = "data/gvcf/{REF}/{sample}.{REF}.g.vcf"
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
        expand("data/processed/scattered_intervals/{REF}/{count}-scattered.interval_list", count = INTERVALS, REF = REF)
    params:
        regions = config.contig_list,
        dir = "data/processed/scattered_intervals/{REF}"
    run:
        shell("gatk SplitIntervals -R {input.ref} -L {params.regions} --scatter-count 200 -O {params.dir}")

# Combine GVCFs with GenomicsDBImport
# sample_map file is samplename with path/to/gvcf in the following line
# scattered interval file created by SpitIntevals

rule combine_gvcfs:
    input:
        gvcfs = expand("data/gvcf/{REF}/{sample}.{REF}.g.vcf", sample = SAMPLE, REF = REF),
        region = "data/processed/scattered_intervals/{REF}/{count}-scattered.interval_list",
        map = config.sample_map
    output:
        directory("data/interm/combined_database_bpres/{REF}/{count}")
    params:
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
        dir = directory("data/interm/combined_database_bpres/{REF}/{count}"),
        ref = config.ref
    output:
        "data/raw/vcf_bpres/{REF}/{count}.raw.{REF}.vcf"
    params:
        db = "gendb://data/interm/combined_database_bpres/{REF}/{count}",
        region = "data/processed/scattered_intervals/{REF}/{count}-scattered.interval_list",
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
