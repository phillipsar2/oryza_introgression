### Liftover vcf coordinates from O. sativa indica to O. glumaepatula

# Prepare maf files downloaded from Ensemble
rule prep_maf:
    input:
        "data/ensemble_alignments/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.maf"
    output:
        swap = temp("data/ensemble_alignments/swap/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.swap.maf")
    run:
        shell("maf-swap {input} > {output.swap}")

rule convert_maf:
    input:
        "data/ensemble_alignments/swap/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.swap.maf"
    output:
        chain = "data/ensemble_alignments/chain/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.chain"
    run:
        shell("grep -v score {input} | maf-convert chain > {output.chain}")

# Combine chain files
rule combine_chain:
    input:
        expand("data/ensemble_alignments/chain/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.{maf}.chain", maf = MAF)
    output:
        "data/ensemble_alignments/merged_chains/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.combined.chain"
    params:
        tmp = "/scratch/aphillip/combine_chain"
    run:
        shell("mkdir -p {params.tmp}")
        shell("chainMergeSort {input} -tempDir={params.tmp} > {output}")
        shell("rm -rf {params.tmp}")

# Subset the vcfs to only the main 12 chromosomes and rename the contigs
# input vcfs need to be bgzipped and tabix indexed
rule rename_vcf:
    input:
        glum = "data/vcf/oryza_glum.vcf.gz",
        indica = "data/vcf/oryza_sativa.half.vcf.gz"
    output:
        glum = "data/vcf/oryza_glum.12chrom.vcf.gz",
        indica = "data/vcf/oryza_sativa.half.12chrom.vcf.gz"
    run:
        shell("bcftools filter {input.glum} -r CM002512.1,CM002513.1,CM002514.1,CM002515.1,CM002516.1,CM002517.1,CM002518.1,CM002519.1,CM002520.1,CM002521.1,CM002522.1,CM002523.1 | \ 
        bcftools annotate --rename-chrs data/vcf/new_glum_names > {output.glum}")
        shell("bcftools filter {input.indica} -r CM000126.1,CM000127.1,CM000128.1,CM000129.1,CM000130.1,CM000131.1,CM000132.1,CM000133.1,CM000134.1,CM000135.1,CM000136.1,CM000137.1 | \
        bcftools annotate --rename-chrs data/vcf/new_indica_names > {output.indica}")

# Liftover the vcf using GATK 4.0
rule liftover:
    input:
        chain = "data/ensemble_alignments/merged_chains/oind_asm465v1.v.oglu_oryza_glumaepatula_v1.5.lastz_net.combined.chain",
        vcf = "data/vcf/oryza_sativa.half.12chrom.vcf",
        ref = "data/genome/GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna"
    output:
        vcf = "data/lifted_vcf/oryza_sativa.half.lifted.vcf",
        reject = "data/lifted_vcf/rejects.vcf"
    run:
        shell("gatk LiftoverVcf \
        -C {input.chain} \
        -I {input.vcf} \
        -O {output.vcf} \
        -R {input.ref} \
        --REJECT {output.reject} \
        --LIFTOVER_MIN_MATCH=0.0 \
        --WRITE_ORIGINAL_POSITION=true \
        --WARN_ON_MISSING_CONTIG=true")
