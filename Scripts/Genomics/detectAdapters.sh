#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

threads=$(( SLURM_CPUS_PER_TASK - 2 ))

module load Python/3.6.4-foss-2018a
module load cutadapt/1.18-foss-2018a-Python-3.6.4

data_name=$1

set -e

Project=b2014286/private

SCRIPTS="/data/pg-molgen/david/scripts"
#TOOLS="/proj/$Project/tools"
FULL_ADAPTERS_FILENAME="candidate-adapters-v4.conf"
DEFAULT_ADAPTERS="/data/pg-molgen/david/scripts/$FULL_ADAPTERS_FILENAME"
if [ -e "./$FULL_ADAPTERS_FILENAME" ] ; then
    DEFAULT_ADAPTERS="./$FULL_ADAPTERS_FILENAME"
fi
OUTPUTDIR="/data/pg-molgen/david/01_postQC-data"
TMP=${SNIC_TMP:-/data/pg-molgen/david/01_postQC-data}
TmpID=${SLURM_JOB_ID:-$$}
# this is available in the module
cutadapt="cutadapt"

Read1=$1"_1.fastq.gz"
Read2=$1"_2.fastq.gz"
ADAPTERS=${3:-$DEFAULT_ADAPTERS}
SUBSET=${4:-1000000}

set +x
if [[ -z "$Read1" || -z "$Read2" || -z "$ADAPTERS" || -z "$TMP" || -z "$SUBSET" ]] ; then
  echo "Usage: $0  read1.fq[.gz] read2.fq[.gz]  [ cutadapt-list-of-adapters ] [ read-subset-size ] "
  echo
  echo "Default file for adapters is $ADAPTERS"
  echo "Default size of read subset is $SUBSET"
  echo
  echo "cutadapt report will be in directory '$OUTPUTDIR', with read1 filename prefix and suffix '.cutReport'"
  echo
  exit 1;
fi
set -x

# /path/to/Reads/Read1.fq.gz
Read1_file=${Read1##*/}                 # Read1.fq.gz
Read1_file_base=${Read1_file%.fastq.gz}                 # Read1.fastq.gz
Read1_file_base=${Read1_file_base%.fq.gz}                 # Read1.fq.gz
Read1_path=${Read1%/*}                  # /path/to/Reads
Read1_dir1=${Read1_path##*/}             # Reads

# /path/to/Reads/Read2.fq.gz
Read2_file=${Read2##*/}                 # Read2.fq.gz
Read2_path=${Read2%/*}                  # /path/to/Reads
Read2_dir1=${Read2_path##*/}             # Reads

set +x
echo
echo ADAPTERS=$ADAPTERS
echo Read1=$Read1
echo Read2=$Read2
echo

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
set -x

mkdir -p $OUTPUTDIR

#  NOTE use of TMP
Read12_interleaved=$TMP/${Read1_file_base}.subset.i.fq.gz

$SCRIPTS/shuffleSequences_fastq.pl --subset $SUBSET $Read1 $Read2 $Read12_interleaved

#  NOTE use of TMP
Cutadapt_output_interleaved=$TMP/${Read1_file_base}.cutadapt.i.fq.gz
Cutadapt_Report=$OUTPUTDIR/$Read1_file_base.cutReport


# exit

echo "cutadapt $Read12_interleaved using adapters in $ADAPTERS ..."
$cutadapt $(<$ADAPTERS) -O 12 -n 2 -o $Cutadapt_output_interleaved $Read12_interleaved > $Cutadapt_Report

if [ -f "$Cutadapt_Report" ] ; then
	rm -f $Cutadapt_output_interleaved $Read12_interleaved
fi

set +x
echo
echo "Final cutadapt output in $Cutadapt_Report"
echo


