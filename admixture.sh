#!/bin/bash

# defaults
prefix=""
threads=1
rename=1

while getopts "s:f:o:t:r" flag; do
  case "${flag}" in
    s) samples="${OPTARG}" ;; # gzipped vcf file of the cohort samples
    f) founders="${OPTARG}" ;; # gzipped vcf file of just the founders
    o) prefix="${OPTARG}" ;; # output prefix
    t) threads="${OPTARG}" ;; # Threads to use
    r) rename=0 ;; # Don't rename chrs before running admix
    *) echo "Invalid option: -${flag}" >&2; exit 1 ;;
  esac
done

# Get the list of sites used for the cohort samples and merge it with the
# Founders to polarize the three source populations
zgrep -v "#" $samples | cut -f 1-2 > sites_list.txt
tabix $samples
tabix $founders
bcftools merge -R sites_list.txt --threads ${threads} \
-Oz -o ${prefix}_merged.vcf.gz $samples $founders

# Use plink to make a bed file, and do some linkage pruning
plink --vcf ${prefix}_merged.vcf.gz --make-bed --allow-extra-chr \
--out ${prefix}_merged_pruned --indep-pairwise 500 10 0.2 --threads ${threads}

# ADMIXTURE requires integer chromosome codes, so replace the character prefix
# in the .bim file
if [ $r -eq 1 ]; then
  sed -i 's/HiC_scaffold_//g' ${prefix}_merged_pruned.bim
fi
admixture ${prefix}_merged_pruned.bed 3 -jX ${threads} > ${prefix}_admix.log

# Plot the proportions using Joanna Meier's R script
Rscript plotADMIXTURE.R
