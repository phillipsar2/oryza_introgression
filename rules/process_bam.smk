rule add_rg:
    input:
        "data/sorted_bam/{sample}.sorted.{REF}.bam"
    output:
        bam = temp(touch("data/interm/addrg/{sample}.rg.{REF}.bam"))
    params:
        tmp = "/scratch/aphillip/addrg/{sample}",
        sample = "{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk AddOrReplaceReadGroups \
        -I={input} \
        -O={output.bam} \
        -RGID=4 \
        -RGLB=lib1 \
        -RGPL=illumina \
        -RGPU=unit1 \
        -RGSM={params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")

rule mark_dups:
    input:
        "data/interm/addrg/{sample}.rg.{REF}.bam"
    output:
        bam = "data/interm/mark_dups/{REF}/{sample}.dedup.{REF}.bam",
        index = "data/interm/mark_dups/{REF}/{sample}.dedup.{REF}.bai",
        metrics = "qc/mark_dup/{sample}_metrics.{REF}.txt"
    params:
        tmp = "/scratch/aphillip/mark_dups/{sample}"
    run:
        # Create a scratch directory
        shell("mkdir -p {params.tmp}")
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={params.tmp}")
        # Remove scratch directory
        shell("rm -rf {params.tmp}")

# Quality metrics with qualimap
rule bamqc:
    input:
        "data/interm/mark_dups/{REF}/{sample}.dedup.{REF}.bam"
    output:
        "reports/bamqc/{sample}_{REF}_stats/qualimapReport.html"
    params:
        dir = "reports/bamqc/{sample}_{REF}_stats"
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 8 \
        -nr 100000 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=64G")
