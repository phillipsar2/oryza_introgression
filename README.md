Evalutating gene flow between Oryza glumaepatula and domesticated rice


The domesticated samples were trimmed with Trimmomatic because they showed high kmer content and low mapping quality in the beginning of the reads.
The glum samples were not trimmed and the qualimap showed good quality alignment.

# Project status

1. SATIVA - The next step was to start joint genotyping across all samples.

2. GLUM - The next step was to count the final # of SNPs and run a PCA with all samples.

After evaluating coverage, I think I need to call SNPs in two sets: all samples (low coverage) and new samples only (7, high coverage). 
These new high coverage samples will need to be downsampled to 1-2X coverage in order to run GL analyses with ANGSD. 

