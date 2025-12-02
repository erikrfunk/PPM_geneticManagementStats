#!/bin/bash

# Extract the individuals from the year desired
input=$1 # path to VCF file: ~/USS/conservation-genetics/PPM/Loco2025/Imputed_vcfs/Cap_PPM_2025_alphapeel_0.98_PL_maf01_snpmiss0.1_nmiss0.3.vcf.gz
samples=$2 # sample list as a single column saved as a .txt

prefix=$(basename $samples .txt)

bcftools view -S $samples $input |\
bcftools view -e MAC==0 -Oz -o ${prefix}.vcf.gz 
