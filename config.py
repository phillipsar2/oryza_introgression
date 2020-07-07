# Genome

## Oryza glumaepatula samples
#ref = "data/genome/GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna"
#ref_int = expand("data/processed/scattered_intervals/{count}-scattered.interval_list", count = INTERVALS)
#int_dir = "data/processed/scattered_intervals"
#int_region = "data/processed/scattered_intervals/{count}-scattered.interval_list"
#haplo_caller = "data/gvcf/{sample}.g.vcf"
#split_int = "data/processed/scattered_intervals/{count}-scattered.interval_list"
#comb_dir = directory("data/interim/combined_database_bpres/{count}")
#db = "gendb://data/interim/combined_database_bpres/{count}"
#joint_out = "data/raw/vcf_bpres/{count}.raw.vcf"
#get_snps_out = "data/raw/vcf_bpres/{count}.raw.snps.vcf"
#hard_filt =  "data/processed/filtered_snps_bpres/{count}.filtered.snps.vcf"
#diag = "reports/filtering/gvcf_{count}.table"
#nocall = "data/processed/filtered_snps_bpres/{count}.filtered.nocall.vcf"
##depth_table = "reports/filtering_bpres/depth_analysis.table"
#filt_dp1 = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.snps.vcf"
#filt_dp2 = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf"
#bgz_vcf = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.snps.vcf.gz"
#combine = "data/processed/filtered_snps_bpres/oryza_glum.vcf.gz"
#whole_dp1 = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.wholegenome.vcf",
#whole_dp2 = "data/processed/filtered_snps_bpres/{count}.filtered.dp3_77.nocall.wholegenome.vcf"
#whole_genome  "data/processed/filtered_snps_bpres/oglum_wholegenome.vcf.gz"

## Oryza sativa samples
ref = "data/genome/GCA_000004655.2_ASM465v1_genomic.fna"
contig_list = "data/genome/indica.contig.list"
sample_map = "data/processed/indica.sample_map"
ref_int = "expand(""data/processed/scattered_intervals/domesticated/{count}-scattered.interval_list"", count = INTERVALS)"
int_region = "data/processed/scattered_intervals/domesticated/{count}-scattered.interval_list" 
int_dir = "data/processed/scattered_intervals/domesticated/trim/"
haplo_caller = "data/gvcf/trim/{sample}.g.vcf"
split_intervals = "data/processed/scattered_intervals/domesticated/{count}-scattered.interval_list"
comb_dir = "data/interim/combined_database_bpres/domesticated/trim/{count}"
db = "gendb://data/interim/combined_database_bpres/domesticated/trim/{count}"
joint_out = "data/raw/vcf_bpres/domesticated/{count}.raw.vcf"
get_snps_out = "data/raw/vcf_bpres/domesticated/trim/{count}.raw.snps.vcf"
hard_filt =  "data/processed/filtered_snps_bpres/domesticated/trim/{count}.filtered.snps.vcf"
diag = "reports/filtering/domesticated/gvcf_{count}.table"
nocall = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.nocall.vcf"
#depth_table = "reports/filtering_bpres/domesticated/depth_analysis.table"
filt_dp1 = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.snps.vcf"
filt_dp2 = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.nocall.snps.vcf"
bgz_vcf = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.nocall.snps.vcf.gz"
combine = "data/processed/filtered_snps_bpres/domesticated/oryza_sativa.half.vcf.gz"
whole_dp1 = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.wholegenome.vcf",
whole_dp2 = "data/processed/filtered_snps_bpres/domesticated/{count}.filtered.dp3_90.nocall.wholegenome.vcf"
whole_genome =  "data/processed/filtered_snps_bpres/domesticated/osativa_wholegenome.vcf.gz"


## Other
addrg_in = "data/mergensort/{sample}.merge.sorted.bam"
mark_in = "data/interm/addrg/{sample}.rg.bam"


