#!/bin/bash

#SBATCH -A b2014286
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 5-00:00:00
#SBATCH -J IndelRealigner_multipleSamples

module load bioinfo-tools
module load bwa/0.7.15
module load java/sun_jdk1.8.0_92
module load GATK/2.7.2
module load picard/1.92
module load samtools/0.1.19




#java -Xmx16g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar    \
#        -T RealignerTargetCreator \
#        -R /proj/b2014286/private/postQC-data/phaseIII/bwa/PombeRef_withAB325691.fa  \
#        $(cat list_bam_to_RealignerTargetCreator)     \
#        -o list_forIndelRealigner.intervals

java -Xmx60g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar    \
     --maxReadsForRealignment 5000000 \
     -T IndelRealigner         \
     -R /proj/b2014286/private/postQC-data/phaseIII/bwa/PombeRef_withAB325691.fa  \
     $(cat list_bam_to_RealignerTargetCreator)   \
    -targetIntervals list_forIndelRealigner.intervals \
     -nWayOut _realigned.bam



