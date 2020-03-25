#!/bin/bash

#SBATCH -A snic2018-8-12
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH -J SLiM


module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4

Nsim=$1
Ngenerations=$2
output=$3

Rscript --vanilla ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_neutral_ABC_time.R $Nsim $Ngenerations $output

