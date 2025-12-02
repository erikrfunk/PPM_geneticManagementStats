#!/bin/bash

input=$1 # vcf.gz
prefix=$(basename $input .vcf.gz)
genomeLength=1812862276 # fixed number used to scale Froh
#zgrep "#CHROM" $1 | cut -f10- | tr "\t" "\n" > temporarySampleList.txt

# Requires the AF info tag to be correct!
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' $input | bgzip -c > freqs.tab.gz
tabix -s1 -b2 -e2 freqs.tab.gz
tabix $input # not sure if this speeds up roh or not...


bcftools roh -o ${prefix}_roh.txt -H 5e-15 \
--threads 12 --AF-file freqs.tab.gz --exclude 'AF<0.05' --exclude 'AF>0.95' \
$input

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
  }' > ${prefix}_individualRohs.txt
