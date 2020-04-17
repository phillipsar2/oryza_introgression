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

# Trim domesticated rice reads with Trimmomatic
rule trim_reads:
    input:
        r1 = "data/rawdata/{sample}_1.fq.gz",
        r2 = "data/rawdata/{sample}_2.fq.gz"
    output:
        f1 ="data/rawdata/trimmed/{sample}.forward.1.fq",
        u1 = "data/rawdata/trimmed/{sample}.unpaired.1.fq",
        r2 = "data/rawdata/trimmed/{sample}.reverse.2.fq",
        u2 = "data/rawdata/trimmed/{sample}.unpaired.2.fq"
    run:
        shell("java -jar /share/apps/Trimmomatic-0.36/trimmomatic.jar PE \
        {input.r1} {input.r2} \
        {output.f1} {output.u1} {output.r2} {output.u2} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

# Align reads to the reference gnome
rule bwa_map:
    input:
        ref = config.ref,
#        r1 = "data/rawdata/{sample}-R1.fastq",
#        r2 = "data/rawdata/{sample}-R2.fastq"
        r1 = "data/rawdata/trimmed/{sample}.forward.1.fq",
        r2 = "data/rawdata/trimmed/{sample}.reverse.2.fq"
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

