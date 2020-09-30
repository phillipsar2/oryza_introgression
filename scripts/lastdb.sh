## Input files
## Script from Asher Hudson

## Prepare both genomes for alignment

module load last

project_path=/group/jrigrp10/aphillip/oryza
ref_file=GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna
ref_name=$( basename $ref_file )
db_prefix=${project_path}/analyses/last/lastdb/${ref_name}-MAM4
ref=${project_path}/data/genome/${ref_file}
lastdb -P0 -uMAM4 -R01 $db_prefix $ref


ref_file2=GCA_000004655.2_ASM465v1_genomic.fna
ref_name2=$( basename $ref_file2 )
db_prefix2=${project_path}/analyses/last/lastdb/${ref_name2}-MAM4
ref2=${project_path}/data/genome/${ref_file2}
lastdb -P0 -uMAM4 -R01 $db_prefix2 $ref2
