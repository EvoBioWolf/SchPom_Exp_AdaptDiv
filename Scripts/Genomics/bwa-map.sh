#!/bin/bash

#SBATCH -A b2010060
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 2:00:00
#SBATCH -J mapping_bwa_cleaning_alignment 

module load bioinfo-tools
module load bwa/0.7.15
module load java/sun_jdk1.8.0_92
module load GATK/2.7.2
module load picard/1.92
module load samtools/0.1.19

Project=b2014286/private
SCRIPTS="/proj/$Project/scripts"
DATADIR="/proj/$Project/postQC-data/phaseIII/cleanreads"
#OUTPUTDIR="/proj/$Project/postQC-data/phaseIII/bwa"
OUTPUTDIR="$PWD"

sample=$1
Read1=$(ls $DATADIR/$sample"_"*_R1_001.postQC.Trimmomatic032.fq.gz)
Read2=$(ls $DATADIR/$sample"_"*_R2_001.postQC.Trimmomatic032.fq.gz)
Orphan=$(ls $DATADIR/$sample"_"*_R1_001.postQC.Trimmomatic032.se.fq.gz)
Reference=$2
Outputfile="$(echo $OUTPUTDIR/$sample.sort.bam)"
ReadGroup="@RG\tID:$1\tSM:Sample\tPL:illumina\tLB:lib1\tPU:unit1" # this read group header will be added to  SAM and BAM and is neccesary for GATK


echo "the files are:"
echo $Read1
echo $Read2
echo $Orphan
echo $Reference
echo $Outputfile

if [[ ! -e $Read1 ]] ; then 
    echo "$Read1 forward reads file NOT FOUND" 
fi
if [[ ! -e $Read2 ]] ; then 
    echo "$Read2 reverse reads file NOT FOUND" 
fi
if [[ ! -e $Orphan ]] ; then 
    echo "$Orphan orphan reads file NOT FOUND" 
fi   
if [[ ! -e $REFDIR/$Reference ]] ; then 
    echo "$Reference reference file  NOT FOUND" 
fi

if [[ -z "$Read1" || -z "$Reference" ]] ; then
  echo "Usage: $0  read1 reference.fa"
  echo
  echo "Directory for output is $OUTPUTDIR"
  echo
  exit 1;
fi
echo "No of reads times 4 Forward"
zcat $Read1 | wc -l
echo "No of reads times 4 Reverse"
zcat $Read2 | wc -l 
echo "No of reads timees 4 Orphan"
zcat $Orphan | wc -l

echo "Reads file fwd:  $Read1"
echo "Reads file rev:  $Read2"
echo "Single end file: $Orphan"
echo "Reference:       $Reference"
echo "Output file:     $Outputfile"

# here there actual alignment of reads to the reference takes place, flag
# -M is needed to mark any multi part alignment as secondary
# then we need to convert the SAM --> BAM and sort the latter
bwa mem -t 3 -R $ReadGroup -M $Reference $Read1 $Read2 | java -jar /sw/apps/bioinfo/picard/1.92/milou/SortSam.jar \
  INPUT=/dev/stdin \
  OUTPUT="$Outputfile" \
  SORT_ORDER=coordinate

bwa mem -t 3 -R $ReadGroup -M $Reference $Orphan  | java -jar /sw/apps/bioinfo/picard/1.92/milou/SortSam.jar \
  INPUT=/dev/stdin \
  OUTPUT="$OUTPUTDIR/$sample.se.sort.bam" \
  SORT_ORDER=coordinate


# the orphan and paired-end mapping BAM files are merged for further analyses into .merg.bam file
samtools merge -@ 3 "$OUTPUTDIR/$sample.merge.bam" "$Outputfile" "$OUTPUTDIR/$sample.se.sort.bam"

rm "$Outputfile" "$OUTPUTDIR/$sample.se.sort.bam"

# create index from bam file
java -jar /sw/apps/bioinfo/picard/1.92/milou/BuildBamIndex.jar \
  I="$OUTPUTDIR/$sample.merge.bam" \
  OUTPUT="$OUTPUTDIR/$sample.merge.bai" \

##  Mark 'artifacts' and correct these (realigner.sh)
###   The deduped bam-file is searched for areas which are likely to have
###   artifacts  from the mapping, using GATK. this generates a '.list' file
java -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar    \
        -T RealignerTargetCreator \
        -R "$Reference"      \
        -I "$OUTPUTDIR/$sample.merge.bam"     \
        -o "$OUTPUTDIR/$sample.target.list"


###   This list will be used to realign the regions around them again with
###   GATK to generate a bam file that has a better alignment than the 
###   initial alignment created with bwa.
java -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar    \
        -T IndelRealigner         \
        -R "$Reference"      \
        -I "$OUTPUTDIR/$sample.merge.bam"     \
        -targetIntervals "$OUTPUTDIR/$sample.target.list" \
        -o "$OUTPUTDIR/$sample.realigned.bam"

rm "$OUTPUTDIR/$sample.merge.bam" "$OUTPUTDIR/$sample.merge.bai"


