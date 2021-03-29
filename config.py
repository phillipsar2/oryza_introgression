### Wildcards ###
#SAMPLE = glob_wildcards("data/samples/{sample}_1.f1.gz").sample

### Sample preparation ###
# Trim Oryza sativa samples

## Align all samples to the O. glumeapatula genome and call SNPs
#ref = "data/genome/GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna"

#contig_list = "data/genome/glum.contig.list"
#sample_map = "data/processed/GLUM.sample_map"


## filtering
#filt_depth = "data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.3dp100.snps.{REF}.vcf"
#dp_nocall = "data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.3dp100.nocall.snps.{REF}.vcf"
#bgzip_vcf = "data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.3dp100.nocall.snps.{REF}.vcf.gz"

## pop structure
chrom_beagle = "data/beagle/{REF}/oryza.{REF}.{chrom}.beagle.gz"
combined_beagle = "data/beagle/{REF}/oryza.{REF}.beagle.gz"
pca = "data/angsd/{REF}/pca/oryza.{REF}.pcangsd.cov"


######################### Sativa ######################
## Align all samples to the O. sativa genome 
ref = "data/genome/GCA_000004655.2_ASM465v1_genomic.fna"

contig_list = "data/genome/indica.contig.list"
sample_map = "data/processed/SAT.sample_map"

## filtering
#filt_depth = "data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.3dp100.snps.{REF}.vcf"
#dp_nocall = "data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.3dp100.nocall.snps.{REF}.vcf"
#bgzip_vcf = "data/processed/filtered_snps_bpres/{REF}/{interval}.filtered.3dp100.nocall.snps.{REF}.vcf.gz"



##### old stuff ####

#combine = "data/processed/filtered_snps_bpres/domesticated/oryza_sativa.half.vcf.gz"
#whole_dp1 = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.wholegenome.vcf",
#whole_dp2 = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.nocall.wholegenome.vcf"
#whole_genome =  "data/processed/filtered_snps_bpres/domesticated/osativa_wholegenome.vcf.gz"



