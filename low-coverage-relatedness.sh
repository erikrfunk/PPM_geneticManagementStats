#!/bin/bash

# conda activate baseclone

# The default will be to calculate likelihoods from all bams,
# but if a vcf with genotype likelihoods already exist for the reference set of
# individuals, provide a vcf file with GL or GP tag for -g,
# Then likelihoods will be calculated for just a target individual
# If providing a vcf, the bam input should be a path to the individual bam file to be tested

# Set defaults
genos="false"
sampleList="false"
out="output"
threads=1

# Get flagged arguments
while getopts "g:s:o:t:" flag; do
  case "${flag}" in
    g) genos="${OPTARG}" ;;
    s) sampleList="${OPTARG}" ;;
    o) out="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    *) echo "Invalid option: -${flag}" >&2; exit 1 ;;
  esac
done

shift "$((OPTIND-1))"

# Get positional arguments
bamlist=$1

# If a list of sample IDs is not provided, get it from the bams
if [[ $sampleList == "false" ]]; then
  if [ -e samplelist.txt ]; then
    echo "samplelist.txt already exists. Delete and rerun or add it with -i."
    exit 1
  fi
  while read -r bam; do
    basename $bam .bam >> samplelist.txt
  done < $bamlist
  sampleList="samplelist.txt"
fi

# Generate the genotype likelihoods
if [[ $genos == false ]]; then
  angsd -b $bamlist -minq 20 -minmapq 20 -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.02 -doVcf 1 -doPost 1 -out ${out} -P $threads
  tabix ${out}.vcf.gz
  bcftools reheader -s $sampleList ${out}.vcf.gz > tmp && mv tmp ${out}.vcf.gz
else
  prefix=$(basename $bamlist .bam)
  # First, generate the target sample vcf using the sites and alleles of the source genos
  zgrep -v "#" $genos | cut -f1,2,4,5 > sites.txt
  angsd sites index sites.txt
  cut -f1 sites.txt | sort -u > regions.txt
  angsd -b $bamlist -sites sites.txt -rf regions.txt -minq 20 -minmapq 20 -gl 2 -domajorminor 3 -domaf 1 -minmaf 0.02 -doVcf 1 -doPost 1 -out ${prefix} -P $threads
  # Then merge the target individual vcf with the source vcf
  # Combine this with the step to convert GLs to PLs and correct PL type to int
  tabix ${prefix}.vcf.gz
  bcftools reheader -s $sampleList ${prefix}.vcf.gz |
  bcftools merge --force-samples -Oz -o ${out}.vcf.gz $genos
fi

# Convert GL to PL using bcftools plugin
# sed command corrects the data type for the PL format tag
tabix ${out}.vcf.gz
bcftools +tag2tag ${out}.vcf.gz -- --gl-to-pl | \
sed 's/ID=PL,Number=G,Type=Float/ID=PL,Number=G,Type=Integer/g' > \
tmp && gzip tmp && mv tmp.gz ${out}_PLs.vcf.gz

# Run ngsRelate using PLs in the VCF file
ngsRelate -h ${out}_PLs.vcf.gz -p $threads -O ${out}_ngsRelate
