#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH -J SLiM_haploid_neutral


module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4

Ngenerations=$1
L=$2
output=$3

echo $Ngenerations $L $output

for sim in {1..10..1}
do
        echo $sim
        Rscript --vanilla ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_h_neu_ABC_time_3MB.R $sim $Ngenerations $L $output.$L"MB"
done


