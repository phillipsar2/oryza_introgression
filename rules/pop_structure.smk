### Rules for population strucutre analysis ###

# convert filtered vcf to beagle format
rule vcf_to_beagle:
    input:
        vcf = "data/processed/filtered_snps_bpres/{REF}/oryza.{REF}.vcf.gz"
    output:
        chrom = config.chrom_beagle
#        all = config.combined_beagle
    params:
        chr = "{chrom}"
    run:
        shell("vcftools --gzvcf {input.vcf} --BEAGLE-PL --stdout --chr {params.chr} > {output.chrom}")
#        shell("cat <(zcat *beagle.txt.gz | head -n1) <(zcat *beagle.txt.gz | grep -v marker) > {output.all}")

# run PCA
rule angsd_pca:
    input:
        config.combined_beagle
    output:
        config.pca
#        pca = "data/angsd_pi/pca/Ppratensis.pcangsd.cov"
    run:
        shell("python tools/pcangsd/pcangsd.py -beagle {input} -o {output}")
