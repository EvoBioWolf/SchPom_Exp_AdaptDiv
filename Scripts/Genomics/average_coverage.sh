#!/bin/bash

#SBATCH -A b2010060
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:20:00
#SBATCH -J average_coverage

module load bioinfo-tools
module load samtools

samtools depth $1.merge_realigned.bam > coverage_values_$1.txt
grep "MT" -w coverage_values_$1.txt | awk '{sum+=$3} END { print "MT_Average:",sum/NR}' > average_coverage_$1.txt
grep "AB325691" -w coverage_values_$1.txt | awk '{sum+=$3} END { print "AB325691_Average:",sum/NR}' >> average_coverage_$1.txt
grep "III" -w coverage_values_$1.txt | awk '{sum+=$3} END { print "III_Average:",sum/NR}' >> average_coverage_$1.txt
grep "II" -w coverage_values_$1.txt | awk '{sum+=$3} END { print "II_Average:",sum/NR}' >> average_coverage_$1.txt
grep "I" -w coverage_values_$1.txt | awk '{sum+=$3} END { print "I_Average:",sum/NR}' >> average_coverage_$1.txt
rm coverage_values_$1.txt




