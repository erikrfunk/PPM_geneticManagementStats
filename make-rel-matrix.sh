#!/bin/bash
input=$1 # vcf file
prefix=$(basename $input .vcf.gz)

plink --vcf ${input} --allow-extra-chr --make-rel
mv plink.rel ${prefix}.rel
mv plink.rel.id ${prefix}.rel.id
sed -i 's/\t/_/g' ${prefix}.rel.id

# Convert the lower tri matrix to pairwise list
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
' ${prefix}.rel.id ${prefix}.rel >> ${prefix}_relatednessList.txt
