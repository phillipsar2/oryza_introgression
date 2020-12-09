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


# Trim Oryza sativa reads
# Trim polyG tails (-g), use 2 threads (-w 2), minimum length is 36 (-l 36), don't filter for qualit (-Q)

#rule fastp_trim:
#    input:
#        r1 = "data/rawdata/{sample}_1.fq.gz",
#        r2 = "data/rawdata/{sample}_2.fq.gz"
#    output:
#        p1 = "data/rawdata/trimmed/{sample}.trim_1.fq.gz",
#        p2 = "data/rawdata/trimmed/{sample}.trim_2.fq.gz",
#    run:
#        shell("fastp -g -w 2 -l 36 -Q \
#        -i {input.r1} -I {input.r2} \
#        -o {output.p1} -O {output.p2} \
#        -h reports/fastp/osativa.html")


# Align all individuals to a reference genome
rule bwa_map:
    input:
        ref = config.ref,
        r1 = "data/samples/final_set/{sample}_1.fq.gz",
        r2 = "data/samples/final_set/{sample}_2.fq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.{REF}bam")
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"

# Takes the input file and stores a sorted version in a different directory.
rule samtools_sort:
    input:
        "data/interm/mapped_bam/{sample}.mapped.{REF}.bam",
    output:
        temp("data/sorted_bam/{sample}.sorted.{REF}.bam"),
    params:
        tmp = "/scratch/aphillip/sort_bam/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} {input} > {output}")
        shell("rm -rf {params.tmp}")

