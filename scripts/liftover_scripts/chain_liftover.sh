### Creating the final chain file and completeing the liftover

module load maven R java GATK

project_path=/group/jrigrp10/aphillip/oryza
ref_file=GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna
ref_name=$( basename $ref_file )
db_prefix=${project_path}/analyses/last/lastdb/${ref_name}-MAM4
ref=${project_path}/data/genome/${ref_file}

query_file_name=GCA_000004655.2_ASM465v1_genomic.fna
query=${project_path}/data/genome/${query_file_name}
query_name=$(basename $query .gz)
query_name=$(basename $query_name .fa)

vcf=${project_path}/data/processed/filtered_snps_bpres/domesticated/oryza_sativa.half.vcf.gz

## Create a chain file with subset of chains that appear in the the net. 
mkdir -p ${project_path}/analyses/last/netchain_subset/
target_net=${project_path}/analyses/last/target_net/"$ref_name"_"$query_name".net
merged_chain=${project_path}/analyses/last/chain/"$ref_name"_"$query_name".all.chain
liftover_chain=${project_path}/analyses/last/netchain_subset/"$ref_name"_"$query_name".liftOver.chain

netChainSubset $target_net $merged_chain  $liftover_chain


## Run liftOver
mkdir -p ${project_path}/analyses/liftover/
lifted_vcf=${project_path}/analyses/liftover/"$ref_name"_"$query_name".vcf
rejects=${project_path}/analyses/liftover/rejects.vcf

gatk LiftoverVcf -I $vcf -C $liftover_chain -O $lifted_vcf -R $ref --REJECT $rejects --LIFTOVER_MIN_MATCH 0.1
