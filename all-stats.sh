#!/bin/bash

# Flagged arguments
samples=false # optional sample list as a single column saved as a .txt to subset the vcf
setPrefix=false # optionally set a prefix different from either input file prefix
genomeLength=1812862276 # fixed number used to scale the het estimate

mkdir outputs # Will hold the final stats files

while getopts "s:p:g:" flag; do
  case "${flag}" in
    s) samples="${OPTARG}" ;;
    p) setPrefix="${OPTARG}" ;;
    g) genomeLength="${OPTARG}" ;;
    *) echo "Invalid option: -${flag}" >&2; exit 1 ;;
  esac
done

shift $((OPTIND - 1))

#Positional Argument
input=$1 # A vcf.gz file - will be subset if providing a sample list: ~/USS/conservation-genetics/PPM/Loco2025/Imputed_vcfs/Cap_PPM_2025_alphapeel_0.98_PL_maf01_snpmiss0.1_nmiss0.3.v$

# Check if vcf needs to be subsampled
if [[ $samples != false ]]; then
  echo "$(date) -- Subsampling the vcf using sample names provided in ${samples}."
  prefix=$(basename $samples .txt)
  bcftools view -S $samples --force-samples $input | bcftools view -e MAC==0 -Oz -o ${prefix}.vcf.gz
  #Reset input to the subsampled vcf
  input=${prefix}.vcf.gz
else
  echo "$(date) -- No sample list provided. Using all samples in the VCF file."
  prefix=$(basename $input .vcf.gz)
fi

# Secondarily check if the prefix needs to be set differently than above
# Note that the name of the subset vcf will stay the same as the sample list name
# so that the vcf file name reflects the samples, and the set prefix reflects the statistics
if [[ $setPrefix != false ]]; then
   echo "$(date) -- Setting file prefix to: ${setPrefix}"
   prefix=${setPrefix}
fi


# Make a samplelist directly from the vcf to avoid issues later on
zgrep "#CHR" ${input} | cut -f10- | tr "\t" "\n" > ${prefix}.ids

# Tally segregating sites as the number of variants in the cohort vcf
echo "$(date) -- Counting the number of segregating sites at the cohort level."
segSites=$(zgrep -v "#" ${input} | wc -l)
echo $segSites > outputs/${prefix}_segregatingSites.txt
echo "$(date) -- Results written to: outputs/${prefix}_segregatingSites.txt"

# A single plink command to calculate multiple stats at once
echo "$(date) -- Running plink to calculate heterozygosity, inbreeding, and relatedness."
plink --vcf ${input} --allow-extra-chr --het --ibc --make-rel --pca --missing

#############################
# Heterozygosity and F
mv plink.het ${prefix}.het
mv plink.ibc ${prefix}.ibc

echo -e "ID\tO(HOM)\tHETrel\tHETwg\tF\tFhat1" > ${prefix}_finalHetTable.txt
awk -v L=$genomeLength -v S=$segSites 'BEGIN{OFS="\t"} FNR>1{print $1"_"$2,$3,($5-$3)/$5,($5-$3)/(L*($5/S)),$6}' ${prefix}.het > ${prefix}.ObsHet
paste ${prefix}.ObsHet <(tail -n +2 ${prefix}.ibc | cut -f 4) >> ${prefix}_finalHetTable.txt
mv ${prefix}_finalHetTable.txt outputs/
echo "$(date) -- Heterozygosity and F written to: outputs/${prefix}_finalHetTable.txt"

#############################
# PCA from plink
# Use the eigenvectors and missingness from plink
Rscript pca-from-plink.R plink.eigenvec plink.imiss outputs/${prefix}

#############################
# Relatedness Matrix
# Awk command converts the lower triangle matrix into a pairwise list
mv plink.rel ${prefix}.rel

echo -e "ID1\tID2\tREL" >> ${prefix}_relatednessList.txt
awk 'NR==FNR{
      a[FNR]=$1;next
    }
    {
      rowID=a[FNR]
      for(i=1; i<=NF; i++) {
        if($i != ""){
        colID=a[i]
        printf "%s\t%s\t%s\n", rowID,colID,$i
        }
      }
    }
' ${prefix}.ids ${prefix}.rel >> ${prefix}_relatednessList.txt
mv ${prefix}_relatednessList.txt outputs/
echo "$(date) -- Relatedness written to: outputs/${prefix}_relatednessList.txt"

# It seems that snpRelate's IBDKING is more in line with expections given
# the pedigree relationships. So run that stat here using the accessory R script
Rscript snpRelate-IBDKING.R ${input} outputs/${prefix}

#############################
# CalculateROH in bcftools
echo "$(date) -- Calculating 1 and 5MB runs of homozygosity using BCFtools."

# Requires the AF info tag to be correct!
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' ${input} | bgzip -c > freqs.tab.gz
tabix -s1 -b2 -e2 freqs.tab.gz
tabix ${input} # not sure if this speeds up roh or not...


bcftools roh -o ${prefix}_roh.txt -H 5e-15 \
--threads 6 --AF-file freqs.tab.gz --exclude 'AF<0.05' --exclude 'AF>0.95' \
${input}

# Convert the single output file into a sample table
grep "RG" ${prefix}_roh.txt | tail -n +2 |\
awk -v L=$genomeLength '
  BEGIN {
    OFS="\t"
  }
  !($2 in one) {one[$2]=0}
  !($2 in five) {five[$2]=0}
  $6>=1000000 {one[$2]+=$6}
  $6>=5000000 {five[$2]+=$6}
  END {
    print "ID","Roh1MB","Froh1MB","Roh5MB","Froh5MB"
    for (sample in one) {
      print sample,one[sample],one[sample]/L,five[sample],five[sample]/L
    }
  }' > outputs/${prefix}_individualRohs.txt
echo "$(date) -- Results written to: outputs/${prefix}_individualRohs.txt"


#############################
echo "$(date) -- Merging all individual statistics into: outputs/${prefix}_allIndStats.txt ..."

awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} {print a[$1],$2,$3,$4,$5}' \
outputs/${prefix}_finalHetTable.txt outputs/${prefix}_individualRohs.txt > outputs/${prefix}_allIndStats.txt

echo "$(date) -- and outputing cohort averages to: outputs/meanCohortStats.txt"
Rscript cohort-means.R outputs/${prefix}_segregatingSites.txt outputs/${prefix}_allIndStats.txt ${prefix}
