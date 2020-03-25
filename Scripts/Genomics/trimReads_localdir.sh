#!/bin/bash

#SBATCH -A b2010060
#SBATCH -J QC_Trimmomatic
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 3:00:00
Threads=16
JavaMem=60g
# Note using 8 threads, half a milou node
# And using 60GB for half a milou node


set -e

module load bioinfo-tools
module load FastQC
module load java

Project=b2014286/private
SCRIPTS="/proj/$Project/scripts"
TOOLS="/proj/$Project/tools"
OUTPUTDIR="$PWD"
TMP=${SNIC_TMP:-$PWD}

trimmomatic="$TOOLS/Trimmomatic-0.32/trimmomatic-0.32.jar"
OperationTag=Trimmomatic032

TmpID=${SLURM_JOB_ID:-$$}

Read1=$1
Read2=$2
AdapterFasta=$3

if [[ -z "$Read1" || -z "$Read2" || -z "$AdapterFasta" ]] ; then
  echo "Usage: $0  read1.fastq[.gz] read2.fastq[.gz]  final-Fasta-of-adapters"
  echo
  echo "Directory for output is $OUTPUTDIR"
  echo
  echo "Adapter file must contain the final set of adapters to be trimmed PLUS THEIR REVERSE-COMPLEMENTS"
  exit 1;
fi
set -x

# /path/to/Reads/Read1.fq.gz
Read1_file=${Read1##*/}                 # Read1.fq.gz
Read1_file_base=${Read1_file%.fastq.gz}                 # Read1.fq.gz
Read1_file_base=${Read1_file_base%.fq.gz}                 # Read1.fq.gz
Read1_path=${Read1%/*}                  # /path/to/Reads
[ "$Read1_path" = "$Read1" ] && Read1_path="."

Read_file_base=${Read1_file_base%.1}  # removes .1 only if present
Read_file_base=${Read1_file_base%_1}  # removes _1 only if present

# /path/to/Reads/Read2.fq.gz
Read2_file=${Read2##*/}                 # Read2.fq.gz
Read2_file_base=${Read2_file%.fastq.gz}                 # Read1.fq.gz
Read2_file_base=${Read2_file_base%.fq.gz}                 # Read1.fq.gz
Read2_path=${Read2%/*}                  # /path/to/Reads
[ "$Read2_path" = "$Read2" ] && Read2_path="."

set +x
if [ ! -e "${Read1}" ] ; then
	if [ -L "${Read1}" ]; then
		echo link to read 1 file ${Read1} is broken
		exit 1
	else
		echo could not find read 1 file ${Read1}
		exit 1
	fi
fi
if [ ! -e "${Read2}" ] ; then
	if [ -L "${Read2}" ]; then
		echo link to read 1 file ${Read2} is broken
		exit 1
	else
		echo could not find read 1 file ${Read2}
		exit 1
	fi
fi

echo
echo AdapterFasta=$AdapterFasta
echo Read1=$Read1
echo Read2=$Read2
echo Base name used for general output: $Read_file_base
echo Base name used for read 1 output: $Read1_file_base
echo Base name used for read 2 output: $Read2_file_base
echo
set -x

#
# FastQC report on raw reads
#

# Editing this to redirect the fastqc results somewhere writable 
RawFastQCResults="$OUTPUTDIR/$Read_file_base/fastqc-results"
mkdir -p $RawFastQCResults

# directories containing the extracted files are ${RawR1%.${suffixIn}}_fastqc and ${RawR2%.${suffixIn}}_fastqc
# within each directory, the file Images/per_base_quality.png holds the read quality plots
RawR1_Qual=$RawFastQCResults/${Read1_file_base}_fastqc/Images/per_base_quality.png
RawR2_Qual=$RawFastQCResults/${Read2_file_base}_fastqc/Images/per_base_quality.png

if [ ! -f "$RawR1_Qual" -o ! -f "$RawR2_Qual" ] ; then
    fastqc --quiet --threads 2 --outdir $RawFastQCResults --extract $Read1 $Read2
else
    echo FastQC already run on both files of raw reads, results in $RawFastQCResults
fi

#
# run trimmomatic
#
mkdir -p $OUTPUTDIR

Trimmomatic_OUT_I_1=$OUTPUTDIR/$Read1_file_base.postQC.$OperationTag.$TmpID.fq.gz
Trimmomatic_OUT_SE_1=$OUTPUTDIR/$Read1_file_base.postQC.$OperationTag.$TmpID.forward.se.fq.gz
Trimmomatic_OUT_I_2=$OUTPUTDIR/$Read2_file_base.postQC.$OperationTag.$TmpID.fq.gz
Trimmomatic_OUT_SE_2=$OUTPUTDIR/$Read2_file_base.postQC.$OperationTag.$TmpID.reverse.se.fq.gz
Trimmomatic_OUT_SE=$OUTPUTDIR/$Read_file_base.postQC.$OperationTag.$TmpID.se.fq.gz

Trimmomatic_OUT_I_1_final=$OUTPUTDIR/$Read1_file_base.postQC.$OperationTag.fq.gz
Trimmomatic_OUT_I_2_final=$OUTPUTDIR/$Read2_file_base.postQC.$OperationTag.fq.gz
Trimmomatic_OUT_SE_final=$OUTPUTDIR/$Read_file_base.postQC.$OperationTag.se.fq.gz

trimmomatic_Output=$OUTPUTDIR/$Read_file_base.postQC.$OperationTag.trimOutput

java -Xmx${JavaMem} -jar $trimmomatic PE -threads $Threads $QualityScoringOption $Read1 $Read2 $Trimmomatic_OUT_I_1 $Trimmomatic_OUT_SE_1 $Trimmomatic_OUT_I_2 $Trimmomatic_OUT_SE_2 ILLUMINACLIP:${AdapterFasta}:1:30:9 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30 2>&1 | tee $trimmomatic_Output 1>&2

zcat $Trimmomatic_OUT_SE_1 $Trimmomatic_OUT_SE_2 | gzip -c - > $Trimmomatic_OUT_SE && rm -f $Trimmomatic_OUT_SE_1 $Trimmomatic_OUT_SE_2

mv -f $Trimmomatic_OUT_I_1 $Trimmomatic_OUT_I_1_final
mv -f $Trimmomatic_OUT_I_2 $Trimmomatic_OUT_I_2_final
mv -f $Trimmomatic_OUT_SE  $Trimmomatic_OUT_SE_final

#
# collect trimmomatic stats
#
Stats=$(grep "^Input Read Pairs:" $trimmomatic_Output)
NInput=$(echo $Stats | cut -d" " -f4)
NPairSurviving=$(echo $Stats | cut -d" " -f7-8)
NForwardSurviving=$(echo $Stats | cut -d" " -f12-13)
NReverseSurviving=$(echo $Stats | cut -d" " -f17-18)
NDropped=$(echo $Stats | cut -d" " -f20-21)

#
# FastQC report on post-QC reads
#

TrimmedFastQCResults="$OUTPUTDIR/fastqc-results"
mkdir -p $TrimmedFastQCResults

fastqc --quiet --threads 3 --outdir $TrimmedFastQCResults --extract $Trimmomatic_OUT_I_1_final $Trimmomatic_OUT_I_2_final $Trimmomatic_OUT_SE_final

# directories containing the extracted files are ${RawR1%.gz}_fastqc and ${RawR2%.gz}_fastqc
# within each directory, the file Images/per_base_quality.png holds the read quality plots
TrimmedR1_Qual=$TrimmedFastQCResults/$Read1_file_base.postQC.$OperationTag.fq_fastqc/Images/per_base_quality.png
TrimmedR2_Qual=$TrimmedFastQCResults/$Read2_file_base.postQC.$OperationTag.fq_fastqc/Images/per_base_quality.png
TrimmedSE_Qual=$TrimmedFastQCResults/$Read_file_base.postQC.$OperationTag.se.fq_fastqc/Images/per_base_quality.png

HTMLResults="$TrimmedFastQCResults/HTML"
mkdir -p $HTMLResults

RawR1_Qual_final=$HTMLResults/$Read1_file_base.raw.per_base_quality.png
RawR2_Qual_final=$HTMLResults/$Read2_file_base.raw.per_base_quality.png
TrimmedR1_Qual_final=$HTMLResults/$Read1_file_base.postQC.$OperationTag.per_base_quality.png
TrimmedR2_Qual_final=$HTMLResults/$Read2_file_base.postQC.$OperationTag.per_base_quality.png
TrimmedSE_Qual_final=$HTMLResults/$Read_file_base.postQC.$OperationTag.se.per_base_quality.png
Final_Quality_Plots=$OUTPUTDIR/$Read_file_base.postQC.$OperationTag.complete.per_base_quality.png
Final_Quality_HTML=$HTMLResults/$Read_file_base.postQC.$OperationTag.complete.html

convert -resize 400 $RawR1_Qual $RawR1_Qual_final
convert -resize 400 $RawR2_Qual $RawR2_Qual_final
convert -resize 400 $TrimmedR1_Qual $TrimmedR1_Qual_final
convert -resize 400 $TrimmedR2_Qual $TrimmedR2_Qual_final
convert -resize 400 $TrimmedSE_Qual $TrimmedSE_Qual_final
convert \( $RawR1_Qual_final $RawR2_Qual_final \
           \( -background white -pointsize 14 label:"$Read_file_base\n\nQuality scoring option: $QualityScoringOption\n\nInput pairs: $NInput\nPairs surviving: $NPairSurviving\nForward surviving: $NForwardSurviving\nReverse surviving: $NReverseSurviving\nDropped: $NDropped\n\n---------\n\nRaw R1           Raw R2\n\nTrimmed R1    Trimmed R2    Trimmed SE" -gravity west \) \
           +append \) \
        \( $TrimmedR1_Qual_final $TrimmedR2_Qual_final $TrimmedSE_Qual_final +append \) \
        -background white -append $Final_Quality_Plots
rm -f $RawR1_Qual_final $RawR2_Qual_final $TrimmedR1_Qual_final $TrimmedR2_Qual_final $TrimmedSE_Qual_final

set +x

cat <<__report__ > $Final_Quality_HTML
<h3>${Read_file_base} FastQC Results</h3>
<img src="$Final_Quality_Plots">
<table border="0" cellspacing="20">
<tr> <th>Sample</th> <th>Input pairs</th> <th>Pairs surviving</th> <th>Forward surviving</th> <th>Reverse surviving</th> <th>Dropped</th> </tr>
<tr> <td>$Read_file_base</td> <td>$NInput</td> <td>$NPairSurviving</td> <td>$NForwardSurviving</td> <td>$NReverseSurviving</td> <td>$NDropped</td> </tr>
</table>
__report__

echo
echo -------------------------
echo
echo Post-QC reads in files:
echo
echo "    $OUTPUTDIR/$Read1_file_base.postQC.$OperationTag.fq.gz"
echo "    $OUTPUTDIR/$Read2_file_base.postQC.$OperationTag.fq.gz"
echo "    $OUTPUTDIR/$Read_file_base.postQC.$OperationTag.se.fq.gz"
echo
echo QC process single-image summary
echo
echo "    $Final_Quality_Plots"
echo
echo QC process HTML summary, do
echo
echo "    firefox $Final_Quality_HTML"
echo
