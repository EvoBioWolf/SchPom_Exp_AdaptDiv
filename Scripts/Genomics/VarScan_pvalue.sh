#!/bin/bash

#SBATCH -A b2010060
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J VarScan

# this script runs VarScan to produce a list of genetic variants from a pooled sample
# example: sbatch /proj/b2014286/private/scripts/VarScan.sh <bam file without extension> <referece fasta file> <min number of reads to call> <min freq> 

module load bioinfo-tools 
module load java
module load samtools 
module load VarScan/2.3.7

Reference=$2
data_folder=$5
bam_file=$data_folder'/'$1.merge_realigned.bam
min_cov=$3
min_fre=$4

echo "bam_file is: " $bam_file
echo "reference is: " $Reference

#samtools mpileup -f $Reference $bam_file | java -jar $VARSCAN_HOME/VarScan.jar pileup2snp --min-coverage $min_cov --min-var-freq $min_fre > snp_$min_cov'_'$min_fre'_'$1.txt
#samtools mpileup -f $Reference $bam_file | java -jar $VARSCAN_HOME/VarScan.jar pileup2indel --min-coverage $min_cov --min-var-freq $min_fre > indel_$min_cov'_'$min_fre'_'$1.txt

samtools mpileup -q 30 -f $Reference $bam_file | java -jar $VARSCAN_HOME/VarScan.jar pileup2snp --p-value 0.05 --vcf-output > snp_pvalue'_'$1.txt
samtools mpileup -q 30 -f $Reference $bam_file | java -jar $VARSCAN_HOME/VarScan.jar pileup2indel --p-value 0.05 --vcf-output > indel_pvalue'_'$1.txt


