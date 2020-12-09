### Step 2: Make alignments
### script from Asher Hudson

#conda install -c mvdbeek ucsc_tools

project_path=/group/jrigrp10/aphillip/oryza
ref_file=GCA_000576495.1_Oryza_glumaepatula_v1.5_genomic.fna
ref_name=$( basename $ref_file )
db_prefix=${project_path}/analyses/last/lastdb/${ref_name}-MAM4
ref=${project_path}/data/genome/${ref_file}

# Which file is the query and which is the ref?
query_file_name=GCA_000004655.2_ASM465v1_genomic.fna
query=${project_path}/data/genome/${query_file_name}
query_name=$(basename $query .gz)
query_name=$(basename $query_name .fa)

mkdir -p ${project_path}/analyses/last/mat/
mat=${project_path}/analyses/last/mat/"$ref_name"_"$query_name".mat

## Find suitable score parameters for aligning the given sequences
last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 $db_prefix $query > $mat

mkdir -p ${project_path}/analyses/last/maf/
maf=${project_path}/analyses/last/maf/"ref_name"_"$query_name".maf

## Align references
lastal -m50 -E0.05 -C2 -p $mat $db_prefix $query > $maf
# m50 makes allowed multiplicity of initial hits 50, each match lengthened
# until it occurs at most this many times
# E is threshold for errors
# C2 makes lastal discard any gapless alignment in 2 or more gapped alignments
# p specifies match/mismatch score matrix
# Then reference db
# And finally the query fasta

mkdir -p ${project_path}/analyses/last/axt
axt=${project_path}/analyses/last/axt/"$ref_name"_"$query_name".axt

## Convert maf to axt
maf-convert axt $maf > $axt

mkdir -p ${project_path}/analyses/last/chain
chain=${project_path}/analyses/last/chain/"$ref_name"_"$query_name".chain
merged_chain=${project_path}/analyses/last/chain/"$ref_name"_"$query_name".all.chain

#conda install -c bioconda ucsc-axtchain

## Chain together axt alignments - this creates the chain
axtChain $axt $ref $query $chain -linearGap=loose -faQ -faT
# axtChain options
# linearGap=loose for distant species
# -faQ The specified qNibDir is a fasta file with multiple sequences for query
# -faT The specified tNibDir is a fasta file with multiple sequences for target

## Combine sorted files into larger sorted file
mkdir -p ${project_path}/analyses/last/chain_merged
merged_chain=${project_path}/analyses/last/chain/"$ref_name"_"$query_name".all.chain
chainMergeSort $chain > $merged_chain

## Print and store total base count in fa files
mkdir -p ${project_path}/analyses/chromsize/
ref_size=${project_path}/analyses/chromsize/"$ref_name".size
query_size=${project_path}/analyses/chromsize/"$query_name".size

if [ ! -f $ref_size ]; then
  faSize $ref -detailed > $ref_size
fi

if [ ! -f $query_size ]; then
  faSize $query -detailed > $query_size
fi
# note: if your chromosomes are id'd with numbers, faSize will add everything
# before .fa in the file name to the beginning of each chromosome id in the size
# file. This causes issues in the next step. Either make all id's non-numeric
# (e.g. 'chr1' instead of '1') or go in to size file and remove prefix in front
# of chromosome ids.

## Remove chains that don't ahve a chance of being netted
mkdir -p ${project_path}/analyses/last/chain_prenet/
chain_prenet=${project_path}/analyses/last/chain_prenet/"$ref_name"_"$query_name".all.pre.chain

chainPreNet $merged_chain $ref_size $query_size $chain_prenet

## Make alignment nets out of chains
mkdir -p ${project_path}/analyses/last/target_net/
mkdir -p ${project_path}/analyses/last/query_net/
target_net=${project_path}/analyses/last/target_net/"$ref_name"_"$query_name".net
query_net=${project_path}/analyses/last/query_net/"$query_name"_"$ref_name".net

chainNet $chain_prenet $ref_size $query_size $target_net $query_net

##### target_net and merged_chain move forward to the next script (chain_liftover.sh)









### ### ### ### ###

## Making 2bit files for netToAxt
mkdir -p ${project_path}/analyses/last/2bit/
query_twobit=${project_path}/analyses/last/2bit/"$query_name".2bit
ref_twobit=${project_path}/analyses/last/2bit/"$ref_name".2bit

if [ ! -f $query_twobit ]; then
  faToTwoBit $query $query_twobit
fi

if [ ! -f $ref_twobit ]; then
  faToTwoBit $ref $ref_twobit
fi

mkdir -p ${project_path}/analyses/last/net_axt/
net_axt=${project_path}/analyses/last/net_axt/"$ref_name"_"$query_name".net.axt

## Convert net (and chain) to axt
netToAxt $target_net $chain_prenet $ref_twobit $query_twobit $net_axt

mkdir -p ${project_path}/analyses/last/net_axt/net_maf
net_maf=${project_path}/analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".net.maf

## convert from axt to maf format
axtToMaf $net_axt $ref_size $query_size $net_maf
