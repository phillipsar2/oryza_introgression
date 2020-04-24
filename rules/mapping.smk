# fastqc to evaluate quality of read data
# Doesn't currently work
#rule fastqc:
#    input:
#        "data/rawdata/{fastq}.fq.gz", 
#    output:
#        directory("reports/fastqc")
#        html="qc/fastqc/{fastq}.html",
#        zip="qc/fastqc/{fastq}_fastqc.zip"
#    params:
#        files = lambda wildcards, input: " ".join(input),
#        dir = "reports/fastqc"
#    shell:
#        "fastqc -o qc/fastqc -f fastq data/rawdata/*.fq.gz"


# Align reads to the reference gnome
rule bwa_map:
    input:
        ref = config.ref,
#        r1 = "data/rawdata/{sample}-R1.fastq",
#        r2 = "data/rawdata/{sample}-R2.fastq"
        r1 = "data/rawdata/{sample}_1.fq.gz",
        r2 = "data/rawdata/{sample}_2.fq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"

# Takes the input file and stores a sorted version in a different directory.
rule samtools_sort:
    input:
        "data/interm/mapped_bam/{sample}.mapped.bam",
    output:
        temp("data/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/scratch/aphillip/sort_bam/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} {input} > {output}")
        shell("rm -rf {params.tmp}")

