#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH -J SLiM_haploid_sel


module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4


Ngenerations=$1
selection=$2
propSused=$3
L=$4
output=$5

echo $Ngenerations $selection $propSused $L $output.$selection.$propSused 

for sim in {1..10..1}
do
        echo $sim
        Rscript --vanilla ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_h_sel_ABC_time_3MB.R $sim $Ngenerations $selection $propSused $L $output.$selection.$propSused.$L"MB" 
done


