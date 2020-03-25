#!/bin/bash

#SBATCH -A b2010060
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -J Filter_repeats_vcf


module load bioinfo-tools vcftools 

input_data=$1
vcf_file1=$1'/'snp_pvalue_$2.txt
vcf_file2=$1'/'indel_pvalue_$2.txt
bed_file=$3

echo "input_data= "$input_data
echo "vcf_file1= "$vcf_file1
echo "vcf_file2= "$vcf_file2
echo "bed_file= "$bed_file

vcftools --vcf $vcf_file1 --out filtered_snp_pvalue_$2 --exclude-bed $bed_file  --removed-sites 
vcftools --vcf $vcf_file2 --out filtered_indel_pvalue_$2 --exclude-bed $bed_file --removed-sites

grep -wvf filtered_snp_pvalue_$2'.removed.sites' $vcf_file1 | sed -e 's/\%//g' > filtered_snp_pvalue_$2.vcf
grep -wvf filtered_indel_pvalue_$2'.removed.sites' $vcf_file2 | sed -e 's/\%//g' > filtered_indel_pvalue_$2.vcf


