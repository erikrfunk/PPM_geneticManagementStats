#!/bin/bash

# defaults
prefix=""
threads=1
rename=1
prune=0

while getopts "s:f:o:t:i:rp" flag; do
  case "${flag}" in
    s) samples="${OPTARG}" ;; # gzipped vcf file of the cohort samples
    f) founders="${OPTARG}" ;; # gzipped vcf file of just the founders
    o) prefix="${OPTARG}" ;; # output prefix
    t) threads="${OPTARG}" ;; # Threads to use
    i) sampmap="${OPTARG}" ;; # Sample map for admixture plot
    r) rename=0 ;; # Don't rename chrs before running admix
    p) prune=1 ;; # Prune for linkage
    *) echo "Invalid option: -${flag}" >&2; exit 1 ;;
  esac
done

# Get the list of sites used for the cohort samples and merge it with the
# Founders to polarize the three source populations
zgrep -v "#" $samples | cut -f 1-2 > sites_list.txt
tabix $samples
tabix $founders
bcftools filter -R sites_list.txt -Oz -o ${founders}_filtered.vcf.gz $founders
tabix ${founders}_filtered.vcf.gz
echo "$(date) -- Merging the founders with the cohort."
bcftools merge --threads ${threads} $samples ${founders}_filtered.vcf.gz | \
bcftools annotate --set-id +'%CHROM\_%POS' -Oz -o ${prefix}_merged.vcf.gz

if [ $prune -eq 1 ]; then
  echo "$(date) -- Pruning sites for linkage."
  plink --vcf ${prefix}_merged.vcf.gz --allow-extra-chr \
  --indep-pairwise 50 5 0.2 --threads ${threads} --out ${prefix}_merged_pruned

  echo "$(date) -- Making plink files."
  plink --vcf ${prefix}_merged.vcf.gz --make-bed \
  --extract ${prefix}_merged_pruned.pruned.in --threads ${threads} \
  --out ${prefix}_merged
else
# Use plink to make a bed file if not pruning
  echo "$(date) -- Making plink files."
  plink --vcf ${prefix}_merged.vcf.gz --make-bed \
  --threads ${threads} --out ${prefix}_merged
fi

# ADMIXTURE requires integer chromosome codes, so replace the character prefix
# in the .bim file
if [ $rename -eq 1 ]; then
  echo "$(date) -- Renaming chromosomes."
  sed -i 's/HiC_scaffold_//g' ${prefix}_merged.bim
fi

echo "$(date) -- Running admixture."
admixture ${prefix}_merged.bed 3 -j${threads} > ${prefix}_admix.log

# Plot the proportions using Joanna Meier's R script
# First order the sample map to match plink's .fam file. Then plot
echo "$(date) -- Sorting the map file."
awk 'NR==FNR{a[$1]=$0;next} {print a[$1"_"$2]}' ${sampmap} <(cut -f1-2 -d" " ${prefix}_merged.fam) > tmp_map.txt
echo "$(date) -- Plotting admixture proportions."
plotADMIXTURE.R -p ${prefix}_merged_plink -i tmp_map.txt -k 3 -m 3 -l DP,SSM,SM,${prefix} -o ${prefix}_admixturePlot

# Clean up
rm tmp_map.txt
mv ${prefix}_merged* outputs/
