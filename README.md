# PPM_geneticManagementStats
Scripts for calculating genetic statistics for the PPM managed population. <br>
1. [All stats](#all-stats)
2. [Admixture](#admixture)
3. [Low coverage relatedness](#low-coverage-relatedness)
4. [Other scripts](#other-scripts)


## All stats
*Usage:*
```
all-stats.sh [options] input.vcf.gz
  -s, sample list
  -p, set a new prefix for output
  -g, set a whole genome length to scale heterozygosity
  -t, threads
```

*Description:*
<br>
This script runs the majority of the stats at once and requires the path to a vcf file.
Three R scripts are required to be in the path or working directory.
These are `cohort-means.R`, `pca-from-plink.R`, and `snpRelate-IBDKING.R`
<br>
<br>
Optionally, an argument can be added that includes a list of samples.
The sample file should be a single column ending in `.txt` and will be used to subset the vcf.
(Subsetting can also be done using the `extract-individuals.sh` script).

Finally, an alternative prefix can be set. Otherwise, the prefix will be determined from either the sample list, or if not subsetting, the vcf prefix.

The `all-stats.sh` script creates a new directory `outputs`, where all final output files get placed. See below for a more detailed
description of each step included in the script.

1. Segregating sites.
<br>
This is calculated as a simple count of the number of variants included in the vcf file.
This assumes the INFO/AC tag is correct in the vcf file.
If the file was subset using a sample list as part of this script, the INFO/AC tag
should have been updated.
This count requires only a single minor allele to be present for a given site to be counted
and is output with the suffix `_segregatingSites.txt`.<br>
- If a stricter filter is preferred, substitute an MAF expression for the MAC expression with your desired minor allele frequency cutoff. For example, a 5% cutoff:<br>
`bcftools view -e MAF<0.05 input.vcf.gz | grep -v "#" | wc -l `

2. Heterozygosity and Inbreeding coefficients: **plink v1.90b4**
<br>
Output is a text file with the suffix `_finalHetTable.txt`.
<br>
Columns are:
- sample name
- observed homozygous sites
- heterozygous sites as total called sites minus observed homozygous sites divided by the hard coded genome length.
- F as calculated by plink with the `-het` flag
- Fhat1 as calculated by plink with the `-ibc` flag

3. Principal Component Analysis: **plink v1.90b4** and R package **ggplot2**
<br>
- Output is a PCA with the suffix `_PCA.pdf`
- The plot uses missingness to color individuals.
- A second PCA is generated that has points labeled with their ID.
- The second PCA uses the suffix `_PCA_labels.pdf`.

4. Relatedness matrix: **plink v1.90b4** and R package **snpRelate vXX**
- Output is a text file with the suffix `_relatednessList.txt`.
- This includes pairwise relatedness coefficients from the `--make-rel` flag.
<br> **This should be validated against known relatives first!*
<br> **The value from snpRelate's IBDKing calculation seems to be more reliable and likely the better result to use.*

5. Runs of homozygosity: **BCFtools v1.9**
- Output is a text file with the suffix `_individualsRohs.txt`.
- This includes for the total tally of roh for 1 and 5MB regions, as well as the froh using the hard coded genome length.
<br> **Note that this makes a potentially large intermediate file. Likely over 1GB per 5 individuals.*

Finally, all the individual-level stats (heterozygosity, F stats, and ROHs) will be merged into a file with the suffix `_allIndStats.txt` and the cohort mean of each statistic into `meanCohortStats.txt`. Note that while all other output files include the cohort prefix in the file name, the cohort means file does not. Instead, each iteration will append the means the to the end of this file so that multiple cohorts can be run within a loop and the file will accumulate the mean stats for each cohort, adding the name of the cohort into the first column.


## Admixture
*Usage:*<br>
```
admixture.sh [options]
  -s, gzipped vcf file of the cohort samples
  -f, gzipped vcf file of the founder samples
  -i, sample map used to plot admixture proportions
  -o, output prefix
  -t, threads
  -r, don't rename chromosomes
  -p, prune for linkage using plink

  Details:
  Both vcf files should be accessible by bcftools. If a file type error gets
  thrown by bcftools, gunzip the file and pipe it into bcftools view adding the
  -Oz flag to let bcftools do the compression. The script will take care of indexing.
```

*Description:*<br>
Along with the cohort vcf created by either the `all-stat.sh` script or the `extract-individuals.sh` script,
use a vcf file of the founder individuals to calculate admixture proportions for the cohort.

A few assumption this script makes: during the generation of the input files using
plink, a linkage pruning is added. The defaults for this step include <br>
 - Testing snps in a window of 50kb
 - Moving in 5 SNP steps
 - Using an r2 threshold of 0.2

 As of now, these settings are hard coded and would need to be changed by altering
 the script. I think it makes sense to add these as arguments to allow flexibility
 and will work on doing that.

Another assumption is about chromosome names. ADMIXTURE won't allow "non-human"
or non-integer chromosome names. The default in this script is to simply remove
the text part of the chromosome name using sed. This is specific to the PPM reference
where chromosomes all begin with "HiC_scaffold_". If a different string is desired
for replacement, that will need to be manually adjusted. If no chromosome renaming
is required, add the -r flag with no argument to switch the behavior off.

Finally, this script will plot the admixture results in R using a script written
by Joanna Meier called `plotADMIXTURE.R`. This is a barplot of admixture
proportions saved as a tiff. You must provide a two-column sample map file with
sample name in column 1 and population name in column 2. If multiple years are
included in the cohort vcf, this could be a map using source population for founders
and year for each cohort. An example that uses the current founders (as of Dec 2025)
is included in the repository as `admixture-samplemap.txt`. The order will be set
by this script to match the .fam file produced by plink. The plink file automatically
splits IDs into Family and ID using an underscore, so this script will paste them
back together using an underscore in order to match the original ID. If this is not
working or not desired let me know.


## Low coverage relatedness
*Usage:*<br>
```
low-coverage-relatedness.sh [options] input_bamlist.txt
  -g, genotype likelihoods VCF file
  -s, sample list
  -o, output prefix
  -t, threads

  Details:
  -g, genotype likelihoods: a VCF file that includes a GL tag (e.g. angsd -doVcf 1)
  -s, sample list: One column text file of the sample names to use in the final file in the same order as the vcf file.
  -o, output: A prefix for the output file
  -t, threads
  ```

*Description:*<br>  
Purpose: To test an individual that might have lost it's p-chip.<br>
Requires: **ANGSD** and **ngsRelate**<br>

When analyzing the reintroduced populations, we are still limited to low coverage
data until imputation can be optimized. The script here will run **ngsRelate** on low
coverage whole genomes using genotype likelihoods. This requires ngsRelate to be
in your path. A version exists at `/home/centos/USS/erik/Programs/ngsRelate/` <br>

NgsRelate can run on a VCF file that includes phred-scaled likelihoods with a PL tag.
By default, this script will generate this file using **ANGSD** and **BCFtools**,
converting the angsd generated GL tag to a PL tag. But because this is a lengthy
process that involves running angsd on all individuals, there is also an option
to use an existing set of files, and simply add a single individual (or more).
Note that if existing files are used, the VCF with genotype likelihoods is a required argument.
ANGSD will automatically rename all individuals to ind#. Providing a sample list
that includes sample names *in the same order* as the vcf file will rename these
before generating the final relatedness table for easier interpretation.


## Other scripts
The other scripts included in this repository are mostly one-off versions of the statistics
included in the `all-stats.sh` script, in case just a single stat is needed. These all slightly differ from the versions implemented in the combined script, and likely need a bit more testing. Let me know if any of these are needed and I can do a quick check.
