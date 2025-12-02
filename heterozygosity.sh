#!/bin/bash

input=$1 # vcf.gz
genomeLength=1812862276 # fixed number used to scale the het estimate (length of 27 chrs)

prefix=$(basename $input .vcf.gz)


plink --vcf ${input} --allow-extra-chr --het --ibc
mv plink.het ${prefix}.het
mv plink.ibc ${prefix}.ibc

echo -e "ID\tO(HOM)\tHET\tF\tFhat1" > ${prefix}_finalHetTable.txt
awk -v L=$genomeLength 'BEGIN{OFS="\t"} FNR>1{print $1"_"$2,$3,($5-$3)/L,$6}' ${prefix}.het > ${prefix}.ObsHet
paste ${prefix}.ObsHet <(tail -n +2 ${prefix}.ibc | cut -f 4) >> ${prefix}_finalHetTable.txt
