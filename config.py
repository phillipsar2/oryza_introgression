# Genome

## Oryza glumaepatula samples
#ref = "data/genome/GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna"
#ref_int = expand("data/processed/scattered_intervals/{count}-scattered.interval_list", count = INTERVALS)
#int_dir = "data/processed/scattered_intervals"
#int_region = "data/processed/scattered_intervals/{count}-scattered.interval_list"
#comb_dir = directory("data/interim/combined_database_bpres/{count}")
#db = "gendb://data/interim/combined_database_bpres/{count}"
#joint_out = "data/raw/vcf_bpres/{count}.raw.vcf"

## Oryza sativa samples
ref = "data/genome/GCA_000004655.2_ASM465v1_genomic.fna"
contig_list = "data/genome/indica.contig.list"
sample_map = "data/processed/indica.sample_map"
ref_int = "expand(""data/processed/scattered_intervals/domesticated/{count}-scattered.interval_list"", count = INTERVALS)"
int_region = "data/processed/scattered_intervals/domesticated/{count}-scattered.interval_list" 
int_dir = "data/processed/scattered_intervals/domesticated/"
comb_dir = "directory(""data/interim/combined_database_bpres/domesticated/{count}"")"
db = "gendb://data/interim/combined_database_bpres/domesticated/{count}"
joint_out = "data/raw/vcf_bpres/domesticated/{count}.raw.vcf"


## Other
addrg_in = "data/mergensort/{sample}.merge.sorted.bam"

mark_in = "data/interm/addrg/{sample}.rg.bam"


