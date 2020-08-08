# ==================================================
# Documents and scripts were written by: Sergio Tusso
# For manuscript: Tusso, et al. (2020). "Gene flow, fitness trade-offs and ancestral variation govern adaptive divergence" (under review). 
# email: situssog@gmail.com
# Jochen Wolf Lab. 
# +++++++++++++++++++++++++++++++++++++++++++++++++

### 
# Raw data for phenotypic and genomic analyses can be found with their corresponding scripts within the folder "Scripts".

### phenotypic analyses:

# all fitness measurements were analyses using the R script fitness_calculation.R
# input data tables in the same folder


#### Genomic analyses:

# Cleaning of raw data
# CLEANUP

detectAdapters_localdir.sh $reads1 $reads2

mv *cutReport ./cutreport/

##  Cleanup of reads using trimmomatic and do quality control of reads

#the next line will produce a fasta file with adaptors found with more than 1000 counts

grep "Trimmed: [0-9][0-9][0-9][0-9]" -B 2 ./cutreport/$sample"_"*.cutReport | cut -d " " -f2,3 | cut -d ";" -f1 | sed 's/^Adapter />/' - | sed "s/'//g" - | sed '/^$/d' | sed '/^-/d' > pombe-adapters_$sample.fa
mv pombe-adapters_00* pombe-Adapters/

# ex. trimReads.sh file_fwd.fq.gz file_rev.fq.gz pombe-adapters.fa

trimReads_localdir.sh $(ls $(echo */Sample_$sample"/*")) ./pombe-Adapters/pombe-adapters_$sample.fa 


# MAPPING
bwa-map.sh $sample PombeRef.fa

# or excluding IndelRealigner:
bwa-map_withoutIndelRealigner.sh $sample PombeRef.fa ; done
printf ' -I %s' *.merge.bam > list_bam_to_RealignerTargetCreator
printf ' -I %s' phaseII/*dedup.bam >> list_bam_to_RealignerTargetCreator
printf ' -I %s' phaseI/bam_withoutDup/*dedup.bam >> list_bam_to_RealignerTargetCreator

java -Xmx8g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar    \
        -T RealignerTargetCreator \
        -R PombeRef_withAB325691.fa  \
        $(cat list_bam_to_RealignerTargetCreator)     \
        -o list_forIndelRealigner.intervals
java -Xmx8g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar    \
        -T IndelRealigner         \
        --maxReadsForRealignment 1000000 \
        -R PombeRef_withAB325691.fa  \
        $(cat list_bam_to_RealignerTargetCreator)   \
        -targetIntervals list_forIndelRealigner.intervals \
        -nWayOut testrealigned.bam

# the previos steps were included in the script: 
IndelRealigner_multipleSamples.sh
sbatch IndelRealigner_multipleSamples.sh list_bam_to_RealignerTargetCreator_realigned.bam





### Measuring mean coverage per chromosome:
sbatch average_coverage.sh $sample 
grep "" average_coverage_000* | sed -e 's/average_coverage_//g' | sed -e 's/\.txt\:/\t/g' | sed -e 's/_Average: /\t/g' > table_coverage_chromosome.txt
sed -i '1s/^/sample\tchromosome\tcoverage\n/' table_coverage_chromosome.txt
#then the table was analysed with the script coverage_plot_script.R

# and for stdev
samtools depth *bamfile*  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
# The denominator needs to be the size of genome where the BAM file was generated against, rather than the number of bases covered at least once (that's what NR is in this context). To get the genome size from the BAM file:
samtools view -H *bamfile* | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'





##### Genotype calling: 
sbatch VarScan_pvalue.sh $sample PombeRef.fa 1 0.01 bamfiles_folder



###### Repeat Mark
RepeatMasker -pa 4 -xsmall -gff -gccalc -s  -species 'schizosaccharomyces pombe' PombeRef.fa > repeatmasker_pombe
# total length of repeats:
cat PombeRef.bed6  | awk '{print $3-$2}' | awk '{sum +=$1} END {print sum}'

awk '(NR>3){$6=$6-1; print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}' PombeRef.fa.out > PombeRef.bed6
sed '1 i\chrom\tchromStart\tchromEnd\trepeat\ttype\tplus' PombeRef.bed6

### Remove variants in repeats
sbatch filter_repeatsFromVCF.sh filderwithVCFfiles $sample PombeRef.bed6 

# The final VCF files were procesed with a R script: GV_comparison_betweenAllsamples.R
# This script produces most of the figures and supplementary figures for genetic data. 
# This includes filtering of GVs, distribution of GVs, measurements of differenciation and diversity, PCAs, size effect plots, statistical test for divergences between ecotypes and simulations and linkage desequilibrium between pairs of variants. 









### Size effect:
# Example with red populations: 
java -Xmx4g -jar snpEff.jar eff -no-downstream -no-upstream Schizosaccharomyces_pombe table_vcf_allGV_formatted_red.vcf > ann_final_table_vcf_allGV_formatted_red.vcf
mv snpEff_genes.txt snpEff_genes_table_vcf_allGV_formatted_red.txt
mv snpEff_summary.html snpEff_summary_table_vcf_allGV_formatted_red.html

grep -v "^#" ann_final_table_vcf_allGV_formatted_red.vcf  | sed 's/|/\t/g' | awk '{print $1,$2,$4,$5,$9,$10,$11 }' | grep -v "CHROM" | sed -e '1iCHROM POS REF ALT Effect Level Gene\' | sed 's/ /\t/g' > ann_final_table_vcf_allGV_formatted_red_ed.txt


 




# SLIM simulations:

# Using the known parameters of pombe I simulated first one single population:
# Ne: 1e6, Window Size: 3Mb, r=1e-6, CloningRate(0.9), u=2e-10

SLiM_1P_neutral_ABC_time.sh num_replicate_populations num_generations outputFile.txt


# Simulations in a 3 MB chromosomes:
SLiM_1P_neutral_ABC_time_3MB.sh num_replicate_populations num_generations outputFile.txt

# Here I simulated 3MB chromosomes but with selection:
# in these simulations I used exponential distribution of benefitial mutations with mean 0.01 (2nd parameter), and 1 every 1000 mutations is benefitial (3rd parameter).
# haploid model with selection:
sbatch SLiM_1P_h_sel_ABC_time_3MB.sh Num_generations mean_sel ratio_sel_neutralMutations 1 replicate_simulation


# haploid model without selection:
# I did two test. with 800 (only phase III) and with 1000 (with Phase II but ignoring that populations were mixed - which reduced frequencies). And just to test I also did 2000
SLiM_1P_h_neu_ABC_time_3MB.sh Num_generations 1 replicate_simulation






