#clean memory
# rm(list=ls(all=TRUE))

#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(Hmisc)

require(gplots) 

# functions:
ggplot_alternative <- function(plot){
  the_data <- data.frame(x <- seq(0, 1, length.out = 100), y = pbeta(x, 1, 10))
  plot
}
nrow0 <- function(x) dim(x)[1]
varname <- function(x) {deparse(substitute(x))}
logit<-function(x){log(x/(1-x))}
logit.reverse<-function(y){exp(y)/(1+exp(y))}
std <- function(x) sd(x)/sqrt(length(x))

# these three are for the PCA section:
library("RColorBrewer")
library(ggfortify)
library(factoextra)


# set working directory
# setwd("C:/Users/sertu336/Desktop/pombe_phase_III/var_pvalue_realigned/filtered_masker")
# setwd("C:/Users/sertu336/Desktop/pombe_phase_III/var_def")
# setwd("C:/Users/sertu336/Desktop/pombe_phase_III/var_set_par/")
setwd("C:/Users/sertu336/LRZ Sync+Share/pombe_phase_III/var_pvalue_realigned/filtered_masker")


# loading coverage data:

mean_coverage_table<-data.frame(read.table("table_coverage_chromosome.txt", T))
total_coverage <- mean_coverage_table %>% 
  group_by(sample,chromosome) %>% spread(chromosome,coverage) %>% 
  mutate(total_mean_coverage=as.integer(((5700000*I)+(4600000*II)+(3500000*III)+(16000*MT))/(13816000)),
         max_coverage=as.integer(total_mean_coverage+(4*sqrt(total_mean_coverage))))

# loading list of VCF files produced with VarScan
list_files_snps<-list.files(path = ".", pattern = "snp_pvalue")
list_files_indels<-list.files(path = ".", pattern = "indel_pvalue")

#### filtering of genetic variances from parental lines (ebc1-4)
#loading files
# natural strains:
samples_excluded<-c("000377", "000378", "000379", "000380", "000381", "000382", "000383", "000384", "000385", "000386", "000387", "000388", "000389", "000394", "000399", "000400")
#single colony samples:
evolved_sc<-c("000395", "000396", "000397", "000398")
# lines with contamination:
lane_5<-c("000228", "000229", "000230", "000231", "000233", "000234", "000251", "000252", "000253", "000254", "000255", "000256", "000257", "000280", "000377", "000378", "000379", "000380", "000381", "000382", "000383", "000384", "000385", "000386")
lane_6<-c("000158", "000237", "000243", "000244", "000272", "000274", "000279", "000305", "000361", "000362", "000387", "000388", "000389", "000390", "000391", "000392", "000393", "000394", "000395", "000396", "000397", "000398", "000399", "000400")
# or including parentals but contaminated:
lane_6<-c("000158", "000237", "000243", "000244", "000272", "000274", "000279", "000305", "000361", "000362", "000387", "000388", "000389", "000394", "000395", "000396", "000397", "000398", "000399", "000400")


exluded_files_snps<-paste("filtered_snp_pvalue_", samples_excluded, ".vcf", sep="")
exluded_files_indels<-paste("filtered_indel_pvalue_", samples_excluded, ".vcf", sep="")
included_files_snps<-list_files_snps[!(list_files_snps %in% exluded_files_snps)]
included_files_indels<-list_files_indels[!(list_files_indels %in% exluded_files_indels)]

read.tables <- function(file.names, ...) {
  require(plyr)
  ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn, T)))
}

table_snps <- read.tables(c(included_files_snps))

table_snps<-separate(table_snps, Filename, c("file1", "file2","file3", "sampletem"), sep = "\\_", remove=FALSE, fill="right") %>%
  separate(sampletem, c("sample", "file4"), sep = "\\.", remove=TRUE, fill="right") %>%
  select(-file1, -file2, -file3, -file4) %>% mutate(sample=as.numeric(sample), 
                                                    type_GV="snps", 
                                                    GV_code=paste(Chrom, "_", Position, "_", VarAllele, sep=""))

table_indels <- read.tables(c(included_files_indels))
table_indels <- separate(table_indels, Filename, c("file1", "file2","file3", "sampletem"), sep = "\\_", remove=FALSE, fill="right") %>%
  separate(sampletem, c("sample", "file4"), sep = "\\.", remove=TRUE, fill="right") %>%
  select(-file1, -file2, -file3, -file4) %>%
  separate(Cons, c("h1", "corrected_VarAllele"), sep = "\\/", remove=FALSE, fill="right") %>% 
  mutate(VarAllele=corrected_VarAllele) %>% 
  mutate(sample=as.numeric(sample), 
         type_GV="indels", 
         GV_code=paste(Chrom, "_", Position, "_", VarAllele, sep="")) %>% 
  select(-h1, -corrected_VarAllele) 


# loading sample information:

sample_information<-data.frame(read.table("phase_III_sample_information_edFractions.txt", T))

sample_information<-merge(sample_information,total_coverage,by="sample", all=TRUE) %>% select(-I,-II,-III, -MT, -AB325691)
table_snps<-merge(sample_information,table_snps,by="sample", all=TRUE) %>% mutate(pop_treatment=paste(population, "_", treatment, sep=""))
table_indels<-merge(sample_information,table_indels,by="sample", all=TRUE) %>% mutate(pop_treatment=paste(population, "_", treatment, sep=""))

table_snps_withoutfiltering<-table_snps
table_indels_withoutfiltering<-table_indels


# table_snps<-table_snps_withoutfiltering
# table_indels<-table_indels_withoutfiltering


#ignore variants with more than 90% or reads supported by one strand
table_snps <- table_snps %>%  filter(c((Reads2Plus/(Reads2Plus+Reads2Minus))>0.1 & (Reads2Plus/(Reads2Plus+Reads2Minus))<0.9) | Chrom=="MT")
table_indels <- table_indels %>%  filter(c((Reads2Plus/(Reads2Plus+Reads2Minus))>0.1 & (Reads2Plus/(Reads2Plus+Reads2Minus))<0.9) | Chrom=="MT")  


### only GV with p-value lower than 0.0001 were included:
table_snps<-table_snps[table_snps$Pvalue<=0.0001,]
table_indels<-table_indels[table_indels$Pvalue<=0.0001,]


### Include only GV with minimum total coverage (15% of mean coverage), and lower than minimun coverage
table_snps<-table_snps %>%  filter((Reads2*100)/VarFreq>total_mean_coverage*0.15)
table_indels<-table_indels %>%  filter((Reads2*100)/VarFreq>total_mean_coverage*0.15)

# Also filtering by max coverage 

table_snps<-table_snps %>%  filter((((Reads2*100)/VarFreq)<(max_coverage*1.3)) | Chrom=="MT") 
table_indels<-table_indels %>%  filter((((Reads2*100)/VarFreq)<(max_coverage*1.3)) | Chrom=="MT")


# ## Coverage plots per treatment
sample_information<-data.frame(read.table("phase_III_sample_information_edFractions.txt", T))
sample_information<-merge(sample_information,total_coverage,by="sample", all=TRUE) %>% select(-I,-II,-III, -MT, -AB325691) 
sample_information_coverage_analyses<- sample_information %>% 
  filter(!(treatment %in% c("Evolved_1", "Evolved_2", "Evolved_3", "Evolved_4")) & 
           time %in% c(-16,-15,-5,-10,-5,0,53,55,57))
time_series_pop <- sample_information %>% filter(time %in% c(0,7)) %>% mutate(pop_treatment=paste(population, "_", treatment, sep="")) %>% 
  select(pop_treatment) %>% unlist()
sample_information_coverage_analyses2<- sample_information %>% mutate(pop_treatment=paste(population, "_", treatment, sep="")) %>% 
  filter(!(treatment %in% c("Evolved_1", "Evolved_2", "Evolved_3", "Evolved_4")) & time>(-1) & 
  pop_treatment %in% time_series_pop)

sample_information_coverage_analyses_fil <- sample_information_coverage_analyses %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))))  


#### Coverage plots:
######

plot<-ggplot(sample_information_coverage_analyses_fil, aes(x=treatment, y=total_mean_coverage)) +
  geom_boxplot() +  
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", binwidth = 4) +
  facet_grid(. ~ time, scale="free",  space = "free") + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"))
plot
ggsave("./plots/01_coverage_meancoverage_all.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)

labels_treatments<-data.frame(treatment=c("P11", "P12", "P13", "P14", "PhaseI_B", "PhaseI_T", "P2P", "P2R", "P3P", "P3R", "Allo_B", "Allo_T", "Para_B", "Para_T", "LM", "SYM", "Evolved_1", "Evolved_2", "Evolved_3", "Evolved_4"), 
                              treatment2=factor(c("Parental 1", "Parental 2", "Parental 5", "Parental 6", "Phase I Bottom", "Phase I Top", "Phase II Parental Pink", "Phase II Parental Red", "Phase III Parental Pink", "Phase III Parental Red", "Allo Bottom", "Allo Top", "Para Bottom", "Para Top", "LM", "Sym", "Evolved_1", "Evolved_2", "Evolved_3", "Evolved_4")))


sample_information_coverage_analyses_fil %>%
  filter(time %in% c(-16,-10,-5,0, 53)) %>% 
  merge(labels_treatments, by="treatment") %>% 
  mutate(treatment2=factor(treatment2, levels=c("Parental 1", "Parental 2", "Parental 5", "Parental 6", "Phase I Bottom", "Phase I Top", "Phase II Parental Pink", "Phase II Parental Red", "Phase III Parental Pink", "Phase III Parental Red", "Allo Bottom", "Allo Top", "Para Bottom", "Para Top", "LM", "Sym", "Evolved_1", "Evolved_2", "Evolved_3", "Evolved_4"))) %>%
  ggplot(aes(x=treatment2, y=total_mean_coverage)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha=0.8) +
  #geom_dotplot(size=2, alpha=0.7, binaxis = "y", stackdir = "center", position = "dodge", binwidth = 4) +
  facet_grid(. ~ time, scale="free",  space = "free") + 
  ylab("Mean Coverage") +
  xlab("") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+ 
  # ggsave("./plots/01_2_coverage_meancoverage_all.png", width = 8, height = 6, dpi = 600)


ggplot(sample_information_coverage_analyses_fil[sample_information_coverage_analyses_fil$time==(-10),], aes(x=population, y=total_mean_coverage)) +
  geom_bar(stat = "identity") +
  #facet_grid(. ~ population, scale="free",  space = "free") + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+ 
  ggsave("./plots/01_coverage_phaseI.png", width = 8, height = 3, dpi = 600)



plot<-ggplot(sample_information_coverage_analyses_fil, aes(x=treatment, y=total_mean_coverage, fill=parental)) +
  geom_boxplot() +  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", binwidth = 4) +
  facet_grid(. ~ time, scale="free",  space = "free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"))
plot
ggsave("./plots/02_coverage_meancoverage_byparental_all.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)

sample_information_coverage_analyses2_fil <- sample_information_coverage_analyses2 %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))))
plot<-ggplot(sample_information_coverage_analyses2_fil, aes(x=as.factor(time), y=total_mean_coverage)) +
  geom_boxplot() +  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", binwidth = 4) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"))
plot
ggsave("./plots/03_coverage_timeSeries.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)

plot<-ggplot(sample_information_coverage_analyses2_fil[sample_information_coverage_analyses2_fil$time!=0,], aes(x=time, y=total_mean_coverage, colour=pop_treatment)) +
  geom_point() + geom_line() + theme(legend.position="none") +
  facet_grid(. ~ pop_treatment) +#, scale="free",  space = "free") 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"))
plot
ggsave("./plots/04_coverage_timeSeries_bypop.png", ggplot_alternative(plot), width = 18, height = 6, dpi = 600)

######

table_total<-rbind(table_snps, table_indels)

#####
# problematic samples due to cross-contamination

table_number<-table_total %>% filter((time>(-1)) | time==(-15)) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1])
plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV)) +
  geom_point()
plot
ggsave("./plots/05_NumberGV_all_samples.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)

plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV)) +
  geom_point() + ylim(250,750)
plot
ggsave("./plots/06_NumberGV_all_samples_1000.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)


# without lains 5 and 6
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) &
                                 time %in% c(0, 53, 55, 57)) %>% select(sample) %>% unique()
table_number<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6))))) %>%
  filter(!(time %in% c(-16,-15,-10,-5,-1))) %>% 
  filter(!(treatment %in% c("Evolved_1", "Evolved_2", "Evolved_3", "Evolved_4", "P11", "P12", "P13", "P14"))) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1], sample=sample[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1])
# write.table(table_number, "samples_withContaminationProblems.txt", quote = F, row.names = F)
plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV)) +
  geom_point() #+ ylim(350,810) + xlim(180,800)
plot
ggsave("./plots/07_NumberGV_withoutLanes5and6.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)


# correlation between number of GV and coverage:
table_number<-table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) &
                                               VarFreq>(-1)) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1]) %>%
  filter(time %in% c(0,53))
plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV, colour=pop_treatment)) +
  geom_point() + theme(legend.position="none") #+ ylim(350,810) + xlim(180,800)
plot
ggsave("./plots/08_NumberGV_allpopulationsT53andT0.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)

# excluding low frequency GVs (lower than 5-6%)
table_number<-table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6))))  &
                                               VarFreq>(6)) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1]) %>%
  filter(time %in% c(0,53))
plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV)) +
  geom_point() + theme(legend.position="none") + geom_smooth(method='lm')
plot
test_correlation<-cor.test(table_number$total_mean_coverage, table_number$number_GV, method = "spearman")
test_correlation
test_correlation$p.value
ggsave("./plots/09_NumberGV_allpopulationsT53_fq_6.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)

plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV, colour=treatment)) +
  geom_point() + theme(legend.position="none") + geom_smooth(method='lm') +
facet_grid(. ~ treatment)#, scale="free",  space = "free")
plot
table_number_treatment <- table_number %>% filter(treatment=="Allo_B")
table_number_treatment <- table_number %>% filter(treatment=="Allo_T")
table_number_treatment <- table_number %>% filter(treatment=="LM")
table_number_treatment <- table_number %>% filter(treatment=="Para_B")
table_number_treatment <- table_number %>% filter(treatment=="Para_T")
table_number_treatment <- table_number %>% filter(treatment=="SYM")
test_correlation<-cor.test(table_number_treatment$total_mean_coverage,table_number_treatment$number_GV, method = "spearman")
test_correlation$p.value
test_correlation
ggsave("./plots/10_NumberGV_allpopulationsT53_fq_6_byTreatments.png", ggplot_alternative(plot), width = 15, height = 6, dpi = 600)

# time series
table_number<-table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) &
                                               pop_treatment %in% as.vector(time_series_pop) &
                                               VarFreq>(-1)) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1]) %>%
  filter(time>(-30))

plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV)) +#, colour=as.factor(time))) +
  geom_point() + geom_smooth(method='lm') #+ theme(legend.position="none")
plot
test_correlation<-cor.test(table_number$total_mean_coverage, table_number$number_GV, method = "spearman")
test_correlation
test_correlation$p.value
ggsave("./plots/11_NumberGV_timeSeries.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)


# time series
table_number<-table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) &
                                               pop_treatment %in% as.vector(time_series_pop) &
                                               VarFreq>(6)) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1])
plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV)) +
  geom_point() + geom_smooth(method='lm') #+ theme(legend.position="none")
plot
test_correlation<-cor.test(table_number$total_mean_coverage, table_number$number_GV, method = "spearman")
test_correlation
test_correlation$p.value
ggsave("./plots/12_NumberGV_timeSeries_fq_6.png", ggplot_alternative(plot), width = 8, height = 6, dpi = 600)


plot<-ggplot(table_number, aes(x=total_mean_coverage, y=number_GV, colour=as.factor(time))) +
  geom_point() + geom_smooth(method='lm') +
  facet_grid(. ~ time)#, scale="free",  space = "free")
plot
table_number_t47 <- table_number %>% filter(time==47)
test_correlation<-cor.test(table_number_t47$total_mean_coverage, table_number_t47$number_GV, method = "spearman")
test_correlation
test_correlation$p.value
ggsave("./plots/13_NumberGV_timeSeries_fq_6_byTime.png", ggplot_alternative(plot), width = 15, height = 6, dpi = 600)


#####
#### comparing parentals in phase I (2500 Vs HiSeq)
table_total %>% filter(time %in% c(-16,-15)) %>% select(pop_treatment) %>% unique
parentals_phaseI<- table_total %>% filter(time %in% c(-16,-15)) %>%
  group_by(pop_treatment, time) %>%
  dplyr::summarize(number_GV = n(), treatment=treatment[1],
                   parental=parental[1], total_mean_coverage=total_mean_coverage[1])
plot<-ggplot(parentals_phaseI, aes(x=as.integer(time), y=number_GV, colour=pop_treatment)) +
  geom_line() #+ ylim(100,2000)
plot

parentals_phaseI<- table_total %>% filter(time %in% c(-16,-15)) %>% group_by(pop_treatment, GV_code) %>%
  dplyr::summarize(number_GV = n())

parentals_filtering2<- table_total %>% filter(time %in% c(-16,-15)) %>% group_by(pop_treatment, GV_code) %>% 
  dplyr::summarize(number_GV = n()) %>% filter(number_GV < 3) %>% ungroup() %>% select(GV_code)
parentals_filtering2

produce_table_pvalue <- function(pv){
  table_snps_filtered<-table_snps_withoutfiltering[table_snps_withoutfiltering$Pvalue<=pv,]
  table_indels_filtered<-table_indels_withoutfiltering[table_indels_withoutfiltering$Pvalue<=pv,]
  table_total_filtered<-rbind(table_snps_filtered, table_indels_filtered)
  parentals_phaseI<- table_total_filtered %>% 
    filter(time %in% c(-16,-15) & 
             GV_code %in% parentals_filtering2$GV_code) %>%
    group_by(pop_treatment, GV_code) %>%
    select(pop_treatment, GV_code, time, VarFreq) %>%
    mutate(time=paste("time_",time,sep="")) %>% group_by(pop_treatment, GV_code, time) %>%
    dplyr::summarize(number_GV = n(), VarFreq_max=max(VarFreq)) %>%
    spread(time,VarFreq_max, fill= NA) %>% #(-1)) %>%
    #filter(`time_-15`>(-1) & `time_-16`>(-1)) %>% 
    gather("time", "VarFreq_max", 4:5) %>% mutate(pvalue_limit=pv)
  return(parentals_phaseI)
}

pvalue<-c(0.01,0.001,0.0001,0.00001,0.000001,0.0000001, 0.00000001)
final_table_parental<-produce_table_pvalue(0.05)
head(final_table_parental)
for (pv in pvalue){
  final_table_parental<-rbind(final_table_parental, produce_table_pvalue(pv))
}

plot<-ggplot(final_table_parental,
             aes(x=VarFreq_max, fill=as.factor(pvalue_limit))) +
  geom_histogram(position="dodge") +
  facet_grid(time ~ pop_treatment, scale="free") + #+ xlim(-2,100)#, scale="free") #,  space = "free") #+ ylim(0,5000)
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"))
plot
ggsave("./plots/18_pvalue_parentals_phaseI.png", ggplot_alternative(plot), width = 15, height = 6, dpi = 600)

plot<-ggplot(final_table_parental[final_table_parental$time=="time_-15",],
             aes(x=VarFreq_max, fill=pop_treatment)) +
  geom_histogram(position="dodge") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(as.factor(pvalue_limit) ~ time) #+ xlim(-2,100)#, scale="free") #,  space = "free") #+ ylim(0,5000)
plot
ggsave("./plots/19_pvalue_HiSeq_phaseI.png", ggplot_alternative(plot), width = 8, height = 8, dpi = 600)

plot<-ggplot(final_table_parental[final_table_parental$time=="time_-16",],
             aes(x=VarFreq_max, fill=pop_treatment)) +
  geom_histogram(position="dodge") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(as.factor(pvalue_limit) ~ time) #+ xlim(-2,100)#, scale="free") #,  space = "free") #+ ylim(0,5000)
plot
ggsave("./plots/20_pvalue_2500_phaseI.png", ggplot_alternative(plot), width = 8, height = 8, dpi = 600)



plot<-ggplot(final_table_parental[final_table_parental$time=="time_-15" &
                                    final_table_parental$pop_treatment=="3_P13",],
             aes(x=VarFreq_max, fill=pop_treatment)) +
  geom_histogram(position="dodge") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(as.factor(pvalue_limit) ~ ., scale="free") #,  space = "free") #+ ylim(0,5000)
plot
ggsave("./plots/15_mutual_GV_parental_phaseI_sample_parentalP13.png", ggplot_alternative(plot), width = 5, height = 8, dpi = 600)

parentals_phaseI<- table_total %>% filter(time %in% c(-16,-15) & GV_code %in% parentals_filtering2$GV_code) %>%
  group_by(pop_treatment, GV_code) %>%
  select(pop_treatment, GV_code, time, VarFreq) %>%
  mutate(time=paste("time_",time,sep="")) %>% group_by(pop_treatment, GV_code, time) %>%
  dplyr::summarize(number_GV = n(), VarFreq_max=max(VarFreq)) %>% spread(time,VarFreq_max, fill=(-1)) %>%
  filter(`time_-15`>(-1) & `time_-16`==(-1))
plot<-ggplot(parentals_phaseI, aes(x=`time_-15`)) +
  geom_histogram() +
  facet_grid(pop_treatment ~ ., scale="free") #,  space = "free") #+ ylim(0,5000)
plot
ggsave("./plots/16_onlyHiSeq_GV_parental_phaseI.png", ggplot_alternative(plot), width = 15, height = 6, dpi = 600)


parentals_phaseI<- table_total %>% filter(time %in% c(-16,-15) & GV_code %in% parentals_filtering2$GV_code) %>%
  group_by(pop_treatment, GV_code) %>%
  select(pop_treatment, GV_code, time, VarFreq) %>%
  mutate(time=paste("time_",time,sep="")) %>% group_by(pop_treatment, GV_code, time) %>%
  dplyr::summarize(number_GV = n(), VarFreq_max=max(VarFreq)) %>% spread(time,VarFreq_max, fill=(-1)) %>%
  filter(`time_-15`==(-1) & `time_-16`>(-1))
plot<-ggplot(parentals_phaseI, aes(x=`time_-16`)) +
  geom_histogram() +
  facet_grid(pop_treatment ~ ., scale="free") #,  space = "free") #+ ylim(0,5000)
plot
ggsave("./plots/17_only2500_GV_parental_phaseI.png", ggplot_alternative(plot), width = 15, height = 6, dpi = 600)


#####
## Parental strains:
treatment_orden <- c("P11","P12","P13","P14","P2P","P3P","P2R","P3R","Allo_T", "Allo_B","Para_T","Para_B","LM","SYM","PhaseI_B","PhaseI_T","Evolved_1","Evolved_2","Evolved_3","Evolved_4")

#table_number<-
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                                       time %in% c(-16) & #, 55, 57) &
                                       VarFreq>(50)) %>%  # & VarFreq<(110)) %>% 
  group_by(GV_code, type_GV, Chrom, Position) %>% 
  dplyr::summarize(Number_parental_strains = n(), 
                   P11_genotype=("P11" %in% treatment)*1,
                   P12_genotype=("P12" %in% treatment)*1, 
                   P13_genotype=("P13" %in% treatment)*1, 
                   P14_genotype=("P14" %in% treatment)*1) %>% 
  ggplot(aes(x=Position, y=Number_parental_strains)) +
  geom_jitter(size=2, alpha=0.3, height = 0.07) + 
  scale_x_continuous(labels = scales::scientific) + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(. ~ Chrom, scale="free") #+
  ggsave("plots/20_1_distribution_GV_parentalStrains.png", width = 8, height = 3, dpi = 400)
  
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16) & #, 55, 57) &
                         VarFreq>(50)) %>%  # & VarFreq<(110)) %>% 
  group_by(GV_code, type_GV, Chrom, Position) %>% 
  dplyr::summarize(Number_parental_strains = n(), 
                   P11_genotype=("P11" %in% treatment)*1,
                   P12_genotype=("P12" %in% treatment)*1, 
                   P13_genotype=("P13" %in% treatment)*1, 
                   P14_genotype=("P14" %in% treatment)*1) %>% 
  write.table(file = "text_files/table_GV_parentalStrains.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16) & #, 55, 57) &
                         VarFreq>(20)) %>%  # & VarFreq<(110)) %>% 
  group_by(GV_code, type_GV, Chrom, Position) %>% 
  dplyr::summarize(Number_parental_strains = n(), 
                   P11_genotype=("P11" %in% treatment)*1,
                   P12_genotype=("P12" %in% treatment)*1, 
                   P13_genotype=("P13" %in% treatment)*1, 
                   P14_genotype=("P14" %in% treatment)*1) %>% 
  write.table(file = "text_files/table_GV_parentalStrains_min20.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#

venn_data_parental<-table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16) & #, 55, 57) &
                         VarFreq>(50)) %>%  # & VarFreq<(110)) %>% 
  group_by(GV_code, type_GV, Chrom, Position) %>% 
  dplyr::summarize(Number_parental_strains = n(), 
                   P1=("P11" %in% treatment)*1,
                   P2=("P12" %in% treatment)*1, 
                   P3=("P13" %in% treatment)*1, 
                   P4=("P14" %in% treatment)*1) 
head(venn_data_parental)
P1<-venn_data_parental %>% filter(P1==1) %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()
P2<-venn_data_parental %>% filter(P2==1) %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()
P3<-venn_data_parental %>% filter(P3==1) %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()
P4<-venn_data_parental %>% filter(P4==1) %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()

venn(list(P1_P=P1,P2_R=P2,P3_P=P3,P4_R=P4))
venn(list(P1=P1,P2=P2))
venn(list(P1=P1,P2=P3))
venn(list(P1=P1,P2=P4))
venn(list(P1=P2,P2=P3))
venn(list(P1=P2,P2=P4))
venn(list(P1=P3,P2=P4))


######### 
# Phase I
# number of genetic variance before correcting for parental variants:


table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(70)) %>% 
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  group_by(population, time) %>% 
  select(population, genotype, time) %>% 
  dplyr::summarize(N_GV = n()) %>% 
  ggplot(aes(x=population, y=N_GV)) +
  geom_bar(stat="identity") + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(. ~ time, scale="free", space="free") #+
  ggsave("plots/20_0_number_GV_PhaseI_minfq70.png", width = 8, height = 3, dpi = 400)

# number of samples per GV
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(70) & 
                         population!="5S3") %>%  
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:26) %>% 
  filter(genotype==1) %>% 
  group_by(GV_code, Chrom, Position) %>% 
  dplyr::summarize(Number_evolved_strains = n()) %>% 
  #filter(Number_evolved_strains>9) %>% data.frame()
  ggplot(aes(x=Position, y=Number_evolved_strains)) +
  geom_jitter(size=2, alpha=0.3, height = 0.07) + 
  scale_x_continuous(labels = scales::scientific) + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(. ~ Chrom, scale="free") #+
  ggsave("plots/20_3_distribution_GV_PhaseI_minfq30.png", width = 8, height = 4, dpi = 400)



# checking if the problem with 5S3 is an deletion:
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) ) %>% 
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  # filter(Position>1316337 & Position<1317995) %>% # ade loci
  select(GV_code, Chrom, Position, P_5S3, P_3) %>% 
  filter(!(P_5S3==0 & P_3==0)) %>%
  mutate(variant=ifelse((P_5S3==0 & P_3==1), ifelse((P_5S3==1 & P_3==1), 0, (-1)), ifelse((P_5S3==1 & P_3==1), 0, (1)))) %>%
  filter(Chrom!="AB325691") %>% 
  ggplot(aes(x=Position, y=variant)) +
  geom_jitter(size=2, alpha=0.3, height = 0.02) + 
  scale_x_continuous(labels = scales::scientific) + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        # legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(. ~ Chrom, scale="free") #+
  ggsave("plots/20_3_distribution_GV_PhaseI_5S3.png", width = 8, height = 4, dpi = 400)

####


table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16) & #, 55, 57) &
                         VarFreq>(50) & 
                         population!="5S3") %>%  
  group_by(GV_code, type_GV, Chrom, Position) %>% 
  dplyr::summarize(Number_parental_strains = n()) %>% 
  #filter(Number_parental_strains %in% c(1,2,3)) %>% 
  ungroup() %>% 
  select(GV_code, Number_parental_strains) -> GV_parental #%>%
  #unlist() %>% as.vector() -> GV_singletons_parental


# Comparing new GV in phase I with calls in parental strains:
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(70) & 
                         population!="5S3") %>%  
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:26) %>% 
  filter(genotype==1) %>% 
  group_by(GV_code, Chrom, Position) %>% 
  dplyr::summarize(Number_evolved_strains = n()) %>% 
  #mutate(singleton_parental=(GV_code %in% GV_singletons_parental)*1) %>% 
  #filter(Number_evolved_strains>9) %>% data.frame()
  filter(Chrom!="AB325691") %>%
  merge(GV_parental, by="GV_code", all.x = TRUE) %>% 
  ggplot(aes(x=Position, y=Number_evolved_strains)) +
  geom_jitter(size=3, alpha=0.3, height = 0.07) + 
  scale_x_continuous(labels = scales::scientific) + 
  scale_y_continuous(breaks=seq(0,15,4), limits=c(0,15)) + 
  expand_limits(x = 0, y = 0) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(Number_parental_strains ~ Chrom, scale="free") #+
  ggsave("plots/20_4_distribution_GV_PhaseI_byparentalVariants_minfq70.png", width = 8, height = 12, dpi = 400)


# New genetic variants by treatment:
treatments_phaseI <- data.frame(population=c("dgP_1B1", "dgP_1B2", "dgP_1B3", "dgP_1S1", "dgP_1S2", "dgP_1S3", "dgP_2B1", "dgP_2B2", "dgP_2B3", "dgP_2S1", "dgP_2S2", "dgP_2S3", "dgP_5B1", "dgP_5B2", "dgP_5B3", "dgP_5S1", "dgP_5S2", "dgP_5S3", "dgP_6B1", "dgP_6B2", "dgP_6B3", "dgP_6S1", "dgP_6S2", "dgP_6S3"), treatment=c("B", "B", "B", "S", "S", "S", "B", "B", "B", "S", "S", "S", "B", "B", "B", "S", "S", "S", "B", "B", "B", "S", "S", "S"))

table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(30) & 
                         population!="5S3") %>%  
  mutate(population=paste0("P_", population), genotype=(VarFreq>30)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:26) %>% 
  merge(treatments_phaseI, by="population") %>% 
  filter(genotype==1) %>% 
  group_by(GV_code, Chrom, Position, treatment) %>% 
  dplyr::summarize(Number_evolved_strains = n()) %>% 
  #mutate(singleton_parental=(GV_code %in% GV_singletons_parental)*1) %>% 
  #filter(Number_evolved_strains>9) %>% data.frame()
  filter(Chrom!="AB325691") %>%
  spread(treatment, Number_evolved_strains, fill=0) %>% 
  ggplot(aes(x=B, y=S)) +
  geom_jitter(size=3, alpha=0.2, height = 0.1, width = 0.1) +
  #scale_x_continuous(labels = scales::scientific) + 
  #scale_y_continuous(breaks=seq(0,15,4), limits=c(0,15)) + 
  #expand_limits(x = 0, y = 0) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  #facet_grid(. ~ Chrom, scale="free") +
  ggsave("plots/20_5_distribution_GV_PhaseI_bytreatment_minfq30.png", width = 4, height = 4, dpi = 400)

table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(70) & 
                         population!="5S3") %>%  
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:26) %>% 
  merge(treatments_phaseI, by="population") %>% 
  filter(genotype==1) %>% 
  group_by(GV_code, Chrom, Position, treatment) %>% 
  dplyr::summarize(Number_evolved_strains = n()) %>% 
  #mutate(singleton_parental=(GV_code %in% GV_singletons_parental)*1) %>% 
  #filter(Number_evolved_strains>9) %>% data.frame()
  filter(Chrom!="AB325691") %>%
  spread(treatment, Number_evolved_strains, fill=0) #%>% 
  write.table(file = "text_files/table_GV_phaseI_bytreatment_minfq70.txt", quote = FALSE, sep = "\t", row.names = FALSE)


# Number of new GVs per population:
table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10)) %>% 
  filter(VarFreq>(70)) %>%
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         #dgP_5S3=(P_5S3!=P_3 & (sum(P_1,P_2,P_3,P_4)==0))*1, 
         dgP_5S3=((P_5S3==1) & ((P_1+P_2+P_3+P_4==0)))*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_5S3, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:27) %>% 
  merge(treatments_phaseI, by="population") %>% 
  filter(genotype==1) %>% 
  group_by(treatment, population) %>% 
  dplyr::summarize(N_newGV = n()) %>% 
  separate(population, c("extra", "population"), sep="_") %>% 
  #summary()
  ggplot(aes(x=population, y=N_newGV)) +
  geom_bar(stat="identity") + 
  ylab("New Genetic variants") +
  xlab("Strain") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  #facet_grid(. ~ Chrom, scale="free") +
  ggsave("plots/20_6_number_new_GV_PhaseI_bystrain_minfq70.png", width = 7, height = 4, dpi = 400)


ancestral_phaseI <- data.frame(population=c("dgP_1B1", "dgP_1B2", "dgP_1B3", "dgP_1S1", "dgP_1S2", "dgP_1S3", "dgP_2B1", "dgP_2B2", "dgP_2B3", "dgP_2S1", "dgP_2S2", "dgP_2S3", "dgP_5B1", "dgP_5B2", "dgP_5B3", "dgP_5S1", "dgP_5S2", "dgP_5S3", "dgP_6B1", "dgP_6B2", "dgP_6B3", "dgP_6S1", "dgP_6S2", "dgP_6S3"), ancestral=c("P3P", "P3P", "P3P", "P3P", "P3P", "P3P", "P3R", "P3R", "P3R", "P3R", "P3R", "P3R", "P3P", "P3P", "P3P", "P3P", "P3P", "P3P", "P3R", "P3R", "P3R", "P3R", "P3R", "P3R"))

table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(70)) %>% 
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         #dgP_5S3=(P_5S3!=P_3 & (sum(P_1,P_2,P_3,P_4)==0))*1, 
         dgP_5S3=((P_5S3==1) & ((P_1+P_2+P_3+P_4==0)))*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_5S3, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:27) %>% 
  merge(ancestral_phaseI, by="population") %>%
  filter(genotype==1) %>% 
  group_by(ancestral, GV_code, Chrom, Position) %>% 
  dplyr::summarize(N_populations = n()) %>% 
  spread(ancestral, N_populations, fill=0) %>% 
  write.table(file = "text_files/table_GV_phaseI_byancestral_minfq70.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         time %in% c(-16,-10) & #, 55, 57) &
                         VarFreq>(70)) %>% 
  mutate(population=paste0("P_", population), genotype=(VarFreq>70)*1) %>% 
  select(GV_code, Chrom, Position, population, genotype) %>% 
  spread(population, genotype, fill=0) %>% 
  mutate(dgP_1B1=(P_1B1!=P_1)*1, 
         dgP_1B2=(P_1B2!=P_1)*1, 
         dgP_1B3=(P_1B3!=P_1)*1, 
         dgP_1S1=(P_1S1!=P_1)*1, 
         dgP_1S2=(P_1S2!=P_1)*1, 
         dgP_1S3=(P_1S3!=P_1)*1, 
         dgP_2B1=(P_2B1!=P_2)*1, 
         dgP_2B2=(P_2B2!=P_2)*1, 
         dgP_2B3=(P_2B3!=P_2)*1, 
         dgP_2S1=(P_2S1!=P_2)*1, 
         dgP_2S2=(P_2S2!=P_2)*1, 
         dgP_2S3=(P_2S3!=P_2)*1, 
         dgP_5B1=(P_5B1!=P_3)*1, 
         dgP_5B2=(P_5B2!=P_3)*1, 
         dgP_5B3=(P_5B3!=P_3)*1, 
         dgP_5S1=(P_5S1!=P_3)*1, 
         dgP_5S2=(P_5S2!=P_3)*1, 
         #dgP_5S3=(P_5S3!=P_3 & (sum(P_1,P_2,P_3,P_4)==0))*1, 
         dgP_5S3=((P_5S3==1) & ((P_1+P_2+P_3+P_4==0)))*1, 
         dgP_6B1=(P_6B1!=P_4)*1, 
         dgP_6B2=(P_6B2!=P_4)*1, 
         dgP_6B3=(P_6B3!=P_4)*1, 
         dgP_6S1=(P_6S1!=P_4)*1, 
         dgP_6S2=(P_6S2!=P_4)*1, 
         dgP_6S3=(P_6S3!=P_4)*1) %>% 
  select(GV_code, Chrom, Position, dgP_1B1, dgP_1B2, dgP_1B3, dgP_1S1, dgP_1S2, dgP_1S3, dgP_2B1, dgP_2B2, dgP_2B3, dgP_2S1, dgP_2S2, dgP_2S3, dgP_5B1, dgP_5B2, dgP_5B3, dgP_5S1, dgP_5S2, dgP_5S3, dgP_6B1, dgP_6B2, dgP_6B3, dgP_6S1, dgP_6S2, dgP_6S3) %>% 
  gather(population, genotype, 4:27) %>% 
  merge(ancestral_phaseI, by="population") %>% 
  filter(genotype==1) %>% 
  group_by(ancestral, GV_code, Chrom, Position) %>% 
  dplyr::summarize(N_populations = n()) %>% 
  filter(ancestral=="P3R" & N_populations>0) %>% dim()


#####



#####
# Phase II

parental_GV<-read.table("text_files/table_GV_parentalStrains.txt", T)
parental_GV_20<-read.table("text_files/table_GV_parentalStrains_min20.txt", T)
phaseI_GV<-read.table("text_files/table_GV_phaseI_byancestral_minfq30.txt", T)

fix_parental_GV_red<-parental_GV %>% filter(P12_genotype==1 & P14_genotype==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()
fix_parental_GV_pink<-parental_GV %>% filter(P11_genotype==1 & P13_genotype==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()

polymorphic_parental_GV_red<-parental_GV %>% filter((P12_genotype+P14_genotype)==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()
polymorphic_parental_GV_pink<-parental_GV %>% filter((P11_genotype+P13_genotype)==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()

#
fix_parental_GV_red_20<-parental_GV_20 %>% filter(P12_genotype==1 & P14_genotype==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()
fix_parental_GV_pink_20<-parental_GV_20 %>% filter(P11_genotype==1 & P13_genotype==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()

polymorphic_parental_GV_red_20<-parental_GV_20 %>% filter((P12_genotype+P14_genotype)==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()
polymorphic_parental_GV_pink_20<-parental_GV_20 %>% filter((P11_genotype+P13_genotype)==1) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()

#

polymorphic_phaseI_red<- phaseI_GV %>% 
  filter(!(GV_code %in% c(fix_parental_GV_red, polymorphic_parental_GV_red_20)) & P3R>0) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()

polymorphic_phaseI_pink<- phaseI_GV %>% 
  filter(!(GV_code %in% c(fix_parental_GV_pink, polymorphic_parental_GV_pink_20)) & P3P>0) %>% 
  ungroup() %>% select(GV_code) %>% unlist() %>% as.vector()

phaseII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           time %in% c(-5) & 
           parental == "P3R") 
phaseII_table_total_red$time_groupGV<-"parental_phaseII"
phaseII_table_total_red$time_groupGV[phaseII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseII_table_total_red$time_groupGV[phaseII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseII_table_total_red$time_groupGV[phaseII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"

phaseII_table_total_red$time_groupGV <- factor(phaseII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII"))

phaseII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           time %in% c(-5) & 
           parental == "P3P") 
phaseII_table_total_pink$time_groupGV<-"parental_phaseII"
phaseII_table_total_pink$time_groupGV[phaseII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseII_table_total_pink$time_groupGV[phaseII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseII_table_total_pink$time_groupGV[phaseII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"

phaseII_table_total_pink$time_groupGV <- factor(phaseII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII"))


phaseII_table_total_red %>% ggplot(aes(x=VarFreq, fill=time_groupGV)) +
  geom_histogram(position="dodge") + #aes(y=..density..), 
  ggtitle("Parental phase II red") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/20_7_distribution_GV_parentalPhaseII_red.png", width = 12, height = 4, dpi = 400)

phaseII_table_total_pink %>% ggplot(aes(x=VarFreq, fill=time_groupGV)) +
  geom_histogram(position="dodge") + #aes(y=..density..), 
  ggtitle("Parental phase II pink") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/20_7_distribution_GV_parentalPhaseII_pink.png", width = 12, height = 4, dpi = 400)                 
                 
head(phaseII_table_total_red)

phaseII_table_total_red %>% group_by(time_groupGV) %>% dplyr::summarize(N_GV = n()) 
phaseII_table_total_pink %>% group_by(time_groupGV) %>% dplyr::summarize(N_GV = n()) 

#####
# phase III parental

polymorphic_parental_phaseII_red <- phaseII_table_total_red %>% 
  filter(time_groupGV=="parental_phaseII") %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()

polymorphic_parental_phaseII_pink <- phaseII_table_total_pink %>% 
  filter(time_groupGV=="parental_phaseII") %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()

phaseIII_parental_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           time %in% c(0) & 
           parental == "P3R") 
phaseIII_parental_table_total_red$time_groupGV<-"parental_phaseIII"
phaseIII_parental_table_total_red$time_groupGV[phaseIII_parental_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_parental_table_total_red$time_groupGV[phaseIII_parental_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_parental_table_total_red$time_groupGV[phaseIII_parental_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_parental_table_total_red$time_groupGV[phaseIII_parental_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_parental_table_total_red$time_groupGV <- factor(phaseIII_parental_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII"))


phaseIII_parental_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           time %in% c(0) & 
           parental == "P3P") 
phaseIII_parental_table_total_pink$time_groupGV<-"parental_phaseIII"
phaseIII_parental_table_total_pink$time_groupGV[phaseIII_parental_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_parental_table_total_pink$time_groupGV[phaseIII_parental_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_parental_table_total_pink$time_groupGV[phaseIII_parental_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_parental_table_total_pink$time_groupGV[phaseIII_parental_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_parental_table_total_pink$time_groupGV <- factor(phaseIII_parental_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII"))



phaseIII_parental_table_total_red %>% ggplot(aes(x=VarFreq, fill=time_groupGV)) +
  geom_histogram(position="dodge") + #aes(y=..density..), 
  ggtitle("Parental phase III red") +
  ylim(c(0,120)) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/20_8_distribution_GV_parentalPhaseIII_red.png", width = 12, height = 4, dpi = 400)

phaseIII_parental_table_total_pink %>% ggplot(aes(x=VarFreq, fill=time_groupGV)) +
  geom_histogram(position="dodge") + #aes(y=..density..), 
  ggtitle("Parental phase III pink") +
  ylim(c(0,120)) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/20_8_distribution_GV_parentalPhaseIII_pink.png", width = 12, height = 4, dpi = 400)


phaseIII_parental_table_total_red %>% group_by(time_groupGV) %>% dplyr::summarize(N_GV = n()) 
phaseIII_parental_table_total_pink %>% group_by(time_groupGV) %>% dplyr::summarize(N_GV = n()) 

phaseIII_parental_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>%
  ggplot(aes(x=VarFreq, fill=time_groupGV)) +
  geom_histogram(position="dodge") + #aes(y=..density..), 
  ggtitle("Parental phase III red") +
  #ylim(c(0,120)) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/20_8_distribution_GV_parentalPhaseIII_red_noFixedParentals.png", width = 12, height = 4, dpi = 400)


phaseIII_parental_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>%
  ggplot(aes(x=VarFreq, fill=time_groupGV)) +
  geom_histogram(position="dodge") + #aes(y=..density..), 
  ggtitle("Parental phase III pink") +
  #ylim(c(0,120)) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/20_8_distribution_GV_parentalPhaseIII_pink_noFixedParentals.png", width = 12, height = 4, dpi = 400)

  
head(phaseIII_parental_table_total_red)
head(phaseIII_parental_table_total_pink)

P3R_var<-phaseIII_parental_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>%
  filter(VarFreq>20) %>%
  ungroup() %>% 
  select(GV_code) %>% 
  unique() %>% 
  unlist() %>% as.vector()

P3P_var<-phaseIII_parental_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>%
  filter(VarFreq>20) %>%
  ungroup() %>% 
  select(GV_code) %>% 
  unique() %>% 
  unlist() %>% as.vector()

venn(list(Red_pop=P3R_var,pink_pop=P3P_var))



#####
# Phase III

#####
# individual plots of distribution of variance per populations in red background
parental_phaseIII_red <- phaseIII_parental_table_total_red %>% 
  filter(time_groupGV=="parental_phaseIII") %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         parental == "P3R" & 
                         time==53) %>% 
  ungroup() %>% 
  select(pop_treatment) %>% 
  unique() %>% unlist() %>% as.vector() -> pop_treatment_list

#pop_treatment_id<-"G8_Para_T"
for (pop_treatment_id in pop_treatment_list){
  print (pop_treatment_id)
  phaseIII_table_total_red<-table_total %>% 
    filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
             (time!=53 | pop_treatment==pop_treatment_id) & 
             (time!=(-16) | VarFreq>70) &
             (time!=(-10) | VarFreq>70) &
             time %in% c(-16,-10,-5,0,53) &
             parental == "P3R" ) 
  phaseIII_table_total_red$time_groupGV<-"phaseIII"
  phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
  phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
  phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
  phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
  phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
  phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))
  phaseIII_table_total_red %>% 
    mutate(time_pop_ID=paste0(time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    spread(time_pop_ID, VarFreq, fill=0) %>% 
    gather("time_pop_ID", "VarFreq", 5:21) %>% 
    separate(time_pop_ID, c("time", "pop_ID"), "__") %>% 
    filter(time==53) %>% 
    # filter(time_groupGV!="fix_parental") %>% # this line is to remove fixed parental variants
    ggplot(aes(x=VarFreq, fill=time_groupGV)) +
    geom_histogram(position="dodge") + #aes(y=..density..), 
    ggtitle(paste0("Phase III red ", pop_treatment_id)) +
    ylim(c(0,130)) +
    theme_classic() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "grey70"), 
          panel.grid.minor = element_line(colour = "grey70"), 
          axis.line = element_line(colour = "black"), 
          #legend.position="bottom", 
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black")) + 
    ggsave(paste0("plots/22_1_red/22_1_distribution_GV_PhaseIII_red_",pop_treatment_id,".png"), width = 10, height = 5, dpi = 400)
    # ggsave(paste0("plots/22_1_red_noParental/22_1_distribution_GV_PhaseIII_red_",pop_treatment_id,".png"), width = 10, height = 5, dpi = 400)  # this line is to remove fixed variants in parentals
}

#####
# individual plots of distribution of variance per populations in pink background
parental_phaseIII_pink <- phaseIII_parental_table_total_pink %>% 
  filter(time_groupGV=="parental_phaseIII") %>% ungroup() %>% 
  select(GV_code) %>% unlist() %>% as.vector()

table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         parental == "P3P" & 
                         time==53) %>% 
  ungroup() %>% 
  select(pop_treatment) %>% 
  unique() %>% unlist() %>% as.vector() -> pop_treatment_list

# pop_treatment_id<-"E1_Allo_T"
for (pop_treatment_id in pop_treatment_list){
  print (pop_treatment_id)
  phaseIII_table_total_pink<-table_total %>% 
    filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
             (time!=53 | pop_treatment==pop_treatment_id) & 
             (time!=(-16) | VarFreq>70) &
             (time!=(-10) | VarFreq>70) &
             time %in% c(-16,-10,-5,0,53) &
             parental == "P3P" ) 
  phaseIII_table_total_pink$time_groupGV<-"phaseIII"
  phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
  phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
  phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
  phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
  phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
  phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))
  phaseIII_table_total_pink %>% 
    mutate(time_pop_ID=paste0(time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    spread(time_pop_ID, VarFreq, fill=0) %>% 
    gather("time_pop_ID", "VarFreq", 5:21) %>% 
    separate(time_pop_ID, c("time", "pop_ID"), "__") %>% 
    filter(time==53) %>% 
    # filter(time_groupGV!="fix_parental") %>% # this line is to remove fixed parental variants
    ggplot(aes(x=VarFreq, fill=time_groupGV)) +
    geom_histogram(position="dodge") + #aes(y=..density..), 
    ggtitle(paste0("Phase III pink ", pop_treatment_id)) +
    ylim(c(0,130)) +
    # ylim(c(0,90)) +
    theme_classic() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "grey70"), 
          panel.grid.minor = element_line(colour = "grey70"), 
          axis.line = element_line(colour = "black"), 
          #legend.position="bottom", 
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black")) + 
    # ggsave(paste0("plots/22_1_pink/22_1_distribution_GV_PhaseIII_pink_",pop_treatment_id,".png"), width = 10, height = 5, dpi = 400)
    ggsave(paste0("plots/22_1_pink_noParental/22_1_distribution_GV_PhaseIII_pink_",pop_treatment_id,".png"), width = 10, height = 5, dpi = 400)  # this line is to remove fixed variants in parentals
}

#####
#### trajectory GV per population: 
# Red:
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         parental == "P3R" & 
                         time==53) %>% 
  ungroup() %>% 
  select(pop_treatment) %>% 
  unique() %>% unlist() %>% as.vector() -> pop_treatment_list

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))


phaseIII_table_total_red %>% 
  filter(time_groupGV=="phaseIII") %>% 
  filter(time %in% c((-16), (-10))) %>% 
  filter(VarFreq>0) %>% 
  ungroup() %>% 
  select(GV_code) %>% 
  unlist() %>% as.vector() -> lowfq_Variants_phaseI_red

# pop_treatment_id<-"G8_Para_T"
# pop_treatment_id<-"B1_Allo_T"

for (pop_treatment_id in pop_treatment_list){
  print (pop_treatment_id)
  phaseIII_table_total_red %>% 
    filter((time!=(-16) | VarFreq>30)) %>% 
    filter((time!=(-10) | VarFreq>30)) %>% 
    mutate(time_pop_ID=paste0(time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    spread(time_pop_ID, VarFreq, fill=0) %>% 
    gather("time_pop_ID", "VarFreq", 5:85) %>% 
    separate(time_pop_ID, c("time", "pop_ID"), "__") %>% 
    filter((time!=(53) | pop_ID==pop_treatment_id)) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=max(VarFreq)) %>% 
    filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>%
    # filter((time_groupGV %in% c("phaseIII"))) %>% 
    mutate(lowFqPhaseI=(GV_code %in% lowfq_Variants_phaseI)*1) %>% 
    # filter(time_groupGV=="phaseI" & time==(-10)) %>% 
    # filter(VarFreq<30)
    ggplot(aes(x=as.factor(as.numeric(time)), y=VarFreq, group=GV_code, colour=time_groupGV)) +
    geom_line() +
    ggtitle(paste0("Red ", pop_treatment_id)) +
    theme_classic() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "grey70"), 
          panel.grid.minor = element_line(colour = "grey70"), 
          axis.line = element_line(colour = "black"), 
          #legend.position="none", 
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black")) +
    facet_grid(time_groupGV ~ .) +
    ggsave(paste0("plots/22_2_trajectory_red/22_2_trajectory_GV_PhaseIII_red_",pop_treatment_id,".png"), width = 10, height = 8, dpi = 400)
  phaseIII_table_total_red %>% 
    mutate(time_pop_ID=paste0(time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    spread(time_pop_ID, VarFreq, fill=0) %>% 
    gather("time_pop_ID", "VarFreq", 5:85) %>% 
    separate(time_pop_ID, c("time", "pop_ID"), "__") %>% 
    filter((time!=(53) | pop_ID==pop_treatment_id)) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=max(VarFreq)) %>% 
    # filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
    filter((time_groupGV %in% c("phaseIII"))) %>% 
    mutate(lowFqPhaseI=(GV_code %in% lowfq_Variants_phaseI_red)*1) %>% 
    # filter(time_groupGV=="phaseI" & time==(-10)) %>% 
    # filter(VarFreq<30)
    ggplot(aes(x=as.factor(as.numeric(time)), y=VarFreq, group=GV_code, colour=time_groupGV)) +
    geom_line() +
    ggtitle(paste0("Red ", pop_treatment_id)) +
    theme_classic() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "grey70"), 
          panel.grid.minor = element_line(colour = "grey70"), 
          axis.line = element_line(colour = "black"), 
          #legend.position="none", 
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black")) +
    facet_grid(time_groupGV ~ lowFqPhaseI) +
    ggsave(paste0("plots/22_2_trajectory_red/lowFq_plots/22_2_trajectory_GV_PhaseIII_red_",pop_treatment_id,"_phaseIIIlowFq.png"), width = 10, height = 8, dpi = 400)
}

# pink populations:
table_total %>% filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
                         parental == "P3P" & 
                         time==53) %>% 
  ungroup() %>% 
  select(pop_treatment) %>% 
  unique() %>% unlist() %>% as.vector() -> pop_treatment_list

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           # (time!=(-16) | VarFreq>30) &
           # (time!=(-10) | VarFreq>30) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))


phaseIII_table_total_pink %>% 
  filter(time_groupGV=="phaseIII") %>% 
  filter(time %in% c((-16), (-10))) %>% 
  filter(VarFreq>0) %>% 
  ungroup() %>% 
  select(GV_code) %>% 
  unlist() %>% as.vector() -> lowfq_Variants_phaseI_pink


for (pop_treatment_id in pop_treatment_list){
  print (pop_treatment_id)
  phaseIII_table_total_pink %>% 
    filter((time!=(-16) | VarFreq>30)) %>% 
    filter((time!=(-10) | VarFreq>30)) %>% 
    mutate(time_pop_ID=paste0(time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    spread(time_pop_ID, VarFreq, fill=0) %>% 
    gather("time_pop_ID", "VarFreq", 5:85) %>% 
    separate(time_pop_ID, c("time", "pop_ID"), "__") %>% 
    filter((time!=(53) | pop_ID==pop_treatment_id)) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=max(VarFreq)) %>% 
    filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>%
    # filter((time_groupGV %in% c("phaseIII"))) %>% 
    mutate(lowFqPhaseI=(GV_code %in% lowfq_Variants_phaseI)*1) %>% 
    # filter(time_groupGV=="phaseI" & time==(-10)) %>% 
    # filter(VarFreq<30)
    ggplot(aes(x=as.factor(as.numeric(time)), y=VarFreq, group=GV_code, colour=time_groupGV)) +
    geom_line() +
    ggtitle(paste0("Pink ", pop_treatment_id)) +
    theme_classic() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "grey70"), 
          panel.grid.minor = element_line(colour = "grey70"), 
          axis.line = element_line(colour = "black"), 
          #legend.position="none", 
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black")) +
    facet_grid(time_groupGV ~ .) +
    ggsave(paste0("plots/22_2_trajectory_pink/22_2_trajectory_GV_PhaseIII_pink_",pop_treatment_id,".png"), width = 10, height = 8, dpi = 400)
  
  phaseIII_table_total_pink %>% 
    mutate(time_pop_ID=paste0(time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    spread(time_pop_ID, VarFreq, fill=0) %>% 
    gather("time_pop_ID", "VarFreq", 5:85) %>% 
    separate(time_pop_ID, c("time", "pop_ID"), "__") %>% 
    filter((time!=(53) | pop_ID==pop_treatment_id)) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time) %>% 
    dplyr::summarize(N_GV=n(), VarFreq=max(VarFreq)) %>% 
    # filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
    filter((time_groupGV %in% c("phaseIII"))) %>% 
    mutate(lowFqPhaseI=(GV_code %in% lowfq_Variants_phaseI_pink)*1) %>% 
    # filter(time_groupGV=="phaseI" & time==(-10)) %>% 
    # filter(VarFreq<30)
    ggplot(aes(x=as.factor(as.numeric(time)), y=VarFreq, group=GV_code, colour=time_groupGV)) +
    geom_line() +
    ggtitle(paste0("Pink ", pop_treatment_id)) +
    theme_classic() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "grey70"), 
          panel.grid.minor = element_line(colour = "grey70"), 
          axis.line = element_line(colour = "black"), 
          #legend.position="none", 
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black")) +
    facet_grid(time_groupGV ~ lowFqPhaseI) +
    ggsave(paste0("plots/22_2_trajectory_pink/lowFq_plots/22_2_trajectory_GV_PhaseIII_pink_",pop_treatment_id,"_phaseIIIlowFq.png"), width = 10, height = 8, dpi = 400)
}

#####
#### Number of GV per time group:
# red

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>% 
  # filter(VarFreq>30) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  group_by(treatment, time_groupGV, pop_ID, time) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  ggplot(aes(pop_ID, N_GV, colour=time_groupGV)) + 
  geom_point(size=3, alpha=0.6) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(. ~ treatment, space="free", scale="free") #+
  # ggsave(paste0("plots/22_3_Number_GV_perGroup_red_minFq30.png"), width = 17, height = 5, dpi = 400)
  ggsave(paste0("plots/22_3_Number_GV_perGroup_red.png"), width = 17, height = 5, dpi = 400)


phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>%  
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>%
  # filter(VarFreq>30) %>%
  filter(time_groupGV!="fix_parental") %>%
  mutate(newGV=factor((time_groupGV=="phaseIII")*1)) %>% 
  group_by(treatment, time_groupGV, pop_ID, time) %>%
  # group_by(treatment, newGV, pop_ID, time) %>%
  dplyr::summarize(N_GV=n()) %>% 
  ggplot(aes(time_groupGV, N_GV, fill=time_groupGV)) +
  # ggplot(aes(newGV, N_GV, fill=newGV)) +
  geom_boxplot() +
  geom_point(size=2, position="jitter", alpha=0.2) +
  scale_y_continuous(breaks=seq(0,150,10)) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("plots/22_3_Number_GV_perGroup_total_red_minFq30.png"), width = 7, height = 5, dpi = 400)
  # ggsave(paste0("plots/22_3_Number_GV_perGroup_total_red.png"), width = 7, height = 5, dpi = 400)
  
  
phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>%   
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(53)) %>% 
  filter(VarFreq>0) %>%
  # filter(VarFreq>20) %>%
  filter(time_groupGV!="fix_parental") %>%
  mutate(newGV=factor(ifelse((time_groupGV=="phaseIII"), "New GV", "Stanting GV"), levels=c("Stanting GV", "New GV"))) %>% 
  # group_by(treatment, time_groupGV, pop_ID, time) %>%
  group_by(newGV, pop_ID, time) %>%
  dplyr::summarize(N_GV=n()) %>% 
  # for the table:
  # group_by(newGV) %>%
  # dplyr::summarize(min_GV=min(N_GV), max_GV=max(N_GV))
  # for the figure:
  ggplot(aes(newGV, N_GV)) +
  # ggplot(aes(newGV, N_GV, fill=newGV)) +
  geom_boxplot() +
  #geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", binwidth = 4) +
  geom_jitter(width = 0.15, alpha=0.8) +
  #geom_point(size=2, position="jitter", alpha=0.6) +
  scale_y_continuous(name="Number Genetic Variants", breaks=seq(0,150,10), limits = c(10, 120)) +
  # scale_y_continuous(name="Number Genetic Variants", breaks=seq(0,150,10), limits = c(0, 45)) +
  facet_grid(time ~ .) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title = element_text(size=14, colour="black"), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
# ggsave(paste0("plots/22_3_2_Number_GV_perGroup_total_red_minFq20.png"), width = 3, height = 5, dpi = 400)
# ggsave(paste0("plots/22_3_2_Number_GV_perGroup_total_red.png"), width = 3, height = 5, dpi = 400)


phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>%   
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(53)) %>% 
  filter(VarFreq>0) %>%
  # filter(VarFreq<20) %>%
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV=="phaseIII") %>%
  # by genetic variants
  # group_by(GV_code) %>% 
  # dplyr::summarize(N_pop=n(), time_groupGV1=time_groupGV[1]) %>% 
  # filter(N_pop==1) %>% 
  # dim()
  group_by(pop_ID) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  summary()
  head()
  filter(N_pop==1) %>% 
  dim()
  
  head()
  

  



# pink
phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

head(phaseIII_table_total_pink)

phaseIII_table_total_pink %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>%
  # filter(VarFreq>30) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  group_by(treatment, time_groupGV, pop_ID, time) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  ggplot(aes(pop_ID, N_GV, colour=time_groupGV)) + 
  geom_point(size=3, alpha=0.6) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(. ~ treatment, space="free", scale="free") #+
  # ggsave(paste0("plots/22_3_Number_GV_perGroup_pink_minFq30.png"), width = 17, height = 5, dpi = 400)
  ggsave(paste0("plots/22_3_Number_GV_perGroup_pink.png"), width = 17, height = 5, dpi = 400)



phaseIII_table_total_pink %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>% 
  # filter(VarFreq>30) %>%
  mutate(newGV=factor((time_groupGV=="phaseIII")*1)) %>% 
  filter(time_groupGV!="fix_parental") %>%
  group_by(treatment, time_groupGV, pop_ID, time) %>%
  # group_by(treatment, newGV, pop_ID, time) %>%
  dplyr::summarize(N_GV=n()) %>% 
  ggplot(aes(time_groupGV, N_GV, fill=time_groupGV)) +
  # ggplot(aes(newGV, N_GV, fill=newGV)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(0,150,10)) +
  geom_point(size=2, position="jitter", alpha=0.2) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("plots/22_3_Number_GV_perGroup_total_pink_minFq30.png"), width = 7, height = 5, dpi = 400)
  # ggsave(paste0("plots/22_3_Number_GV_perGroup_total_pink.png"), width = 7, height = 5, dpi = 400)

phaseIII_table_total_pink %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>%   
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0, 53)) %>% 
  # filter(VarFreq>0) %>%
  filter(VarFreq>20) %>%
  filter(time_groupGV!="fix_parental") %>%
  mutate(newGV=factor(ifelse((time_groupGV=="phaseIII"), "New GV", "Stanting GV"), levels=c("Stanting GV", "New GV"))) %>% 
  # group_by(treatment, time_groupGV, pop_ID, time) %>%
  group_by(newGV, pop_ID, time) %>%
  dplyr::summarize(N_GV=n()) %>% 
  # for the table:
  # group_by(newGV) %>%
  # dplyr::summarize(min_GV=min(N_GV), max_GV=max(N_GV))
  # for the figure:
  ggplot(aes(newGV, N_GV)) +
  # ggplot(aes(newGV, N_GV, fill=newGV)) +
  geom_boxplot() +
  #geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", binwidth = 4) +
  geom_jitter(width = 0.15, alpha=0.8) +
  #geom_point(size=2, position="jitter", alpha=0.6) +
  # scale_y_continuous(name="Number Genetic Variants", breaks=seq(0,150,10), limits = c(10, 120)) +
  scale_y_continuous(name="Number Genetic Variants", breaks=seq(0,150,10), limits = c(0, 45)) +
  facet_grid(time ~ .) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title = element_text(size=14, colour="black"), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("plots/22_3_2_Number_GV_perGroup_total_pink_minFq20.png"), width = 3, height = 5, dpi = 400)
  # ggsave(paste0("plots/22_3_2_Number_GV_perGroup_total_pink.png"), width = 3, height = 5, dpi = 400)


#####
# Number of GVs per treatment - boxplot:
# red

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>% 
  # filter(VarFreq>30) %>%
  filter(time_groupGV!="fix_parental") %>% 
  group_by(treatment, pop_ID, time) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  # filter(N_GV<30)
  # ggplot(aes(N_GV)) +
  # geom_histogram()
  ggplot(aes(treatment, N_GV)) + 
  geom_boxplot() +
  ggtitle("22_3_Number_GV_pertreatment_total_red") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("plots/22_3_Number_GV_pertreatment_total_red.png"), width = 8, height = 7, dpi = 400)


phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>%
  # filter(VarFreq>30) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  group_by(treatment, time_groupGV, pop_ID, time) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  ggplot(aes(treatment, N_GV, fill=time_groupGV)) + 
  geom_boxplot() +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("plots/22_3_Number_GV_pertreatment_red.png"), width = 8, height = 7, dpi = 400)
  # ggsave(paste0("plots/22_3_Number_GV_pertreatment_red_minFq30.png"), width = 8, height = 7, dpi = 400)



# Pink
phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_pink %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>% 
  # filter(VarFreq>30) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  group_by(treatment, pop_ID, time) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  # filter(N_GV<25)
  # ggplot(aes(N_GV)) +
  # geom_histogram()
  ggplot(aes(treatment, N_GV)) + 
  geom_boxplot() +
  ggtitle("22_3_Number_GV_pertreatment_total_pink") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("plots/22_3_Number_GV_pertreatment_total_pink.png"), width = 8, height = 7, dpi = 400)

phaseIII_table_total_pink %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  filter(VarFreq>0) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  group_by(treatment, time_groupGV, pop_ID, time) %>% 
  dplyr::summarize(N_GV=n()) %>% 
  ggplot(aes(treatment, N_GV, fill=time_groupGV)) + 
  geom_boxplot() +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("plots/22_3_Number_GV_pertreatment_pink.png"), width = 8, height = 7, dpi = 400)

#

#####
# Genetic diversity:

# red populations:
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53, 55, 57) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>%
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(0, 53, 55, 57)) %>% 
  filter(time %in% c(53, 55, 57)) %>% 
  # filter(VarFreq>5) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:112) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(VarFreq!=0) %>% head()
  # mutate(pi=(VarFreq*(100-VarFreq))/(((100-1)*100)/2)) %>%
  mutate(pi=2*(VarFreq/100)*((100-VarFreq)/100)) %>%
  # mutate(pi=(((VarFreq-1)*VarFreq)/2)/(((100-1)*100)/2)) %>% 
  # filter(time_groupGV=="phaseIII") %>%
  group_by(treatment, time, pop_ID) %>% 
  dplyr::summarize(N_GV=n(), total_pi=mean(pi)) %>% 
  ungroup() %>% 
  mutate(treatment=factor(treatment, 
                          levels=c("P3R", "Allo_B", "Allo_T", "Para_B", "Para_T", "LM", "SYM"))) %>%
  mutate(fraction=if_else(time %in% c(0,53), "Pool", if_else(time==55, "Top", "Bottom"))) %>%
  ggplot(aes(treatment, total_pi, fill=fraction)) + 
  # geom_boxplot(fill="gray70") +
  geom_boxplot() +
  # geom_point() +
  # geom_point(position="jitter", alpha=0.7, aes(colour=treatment)) + #"black") +
  ggtitle("Pi per population - red") +
  scale_fill_manual(values=c("red3", "orange", "cornflowerblue")) +
  # ggtitle("Pi per population - only new GVs phase III - red") +
  # ylim(c(0.015,0.045)) +
  # ylim(c(0,0.02)) + # this one is for new GVs
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) + 
# ggsave(paste0("plots/22_3_pi_pertreatment_total_red.png"), width = 8, height = 7, dpi = 400)
  ggsave(paste0("plots/22_3_pi_pertreatment_newGV_red.png"), width = 8, height = 7, dpi = 400)


test_table_pi<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>%
  filter(time %in% c(0, 53)) %>% 
  # filter(VarFreq>15) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  mutate(pi=(VarFreq*(100-VarFreq))/(((100-1)*100)/2)) %>% 
  # filter(time_groupGV=="phaseIII") %>%
  group_by(treatment, time, pop_ID) %>% 
  dplyr::summarize(N_GV=n(), total_pi=mean(pi)) %>% 
  ungroup() %>% 
  mutate(treatment=factor(treatment, 
                          levels=c("Allo_B", "Allo_T", "Para_B", "Para_T", "LM", "SYM", "P3R"))) %>%
  filter(treatment!="P3R") %>% 
  data.frame()

head(test_table_pi)
model<-lm(total_pi~treatment, test_table_pi) # to be run with phylogenetic control
summary(model)


# pink populations:
phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53, 55, 57) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>%
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(0, 53, 55, 57)) %>% 
  filter(time %in% c(53, 55, 57)) %>% 
  # filter(VarFreq>5) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:112) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(VarFreq!=0) %>% head()
  # mutate(pi=(VarFreq*(100-VarFreq))/(((100-1)*100)/2)) %>%
  mutate(pi=2*(VarFreq/100)*((100-VarFreq)/100)) %>%
  # mutate(pi=(((VarFreq-1)*VarFreq)/2)/(((100-1)*100)/2)) %>% 
  # filter(time_groupGV=="phaseIII") %>%
  group_by(treatment, time, pop_ID) %>% 
  dplyr::summarize(N_GV=n(), total_pi=mean(pi)) %>% 
  ungroup() %>% 
  mutate(treatment=factor(treatment, 
                          levels=c("P3P", "Allo_B", "Allo_T", "Para_B", "Para_T", "LM", "SYM"))) %>%
  mutate(fraction=if_else(time %in% c(0,53), "Pool", if_else(time==55, "Top", "Bottom"))) %>%
  ggplot(aes(treatment, total_pi, fill=fraction)) + 
  # geom_boxplot(fill="gray70") +
  geom_boxplot() +
  # geom_point() +
  # geom_point(position="jitter", alpha=0.7, aes(colour=treatment)) + #"black") +
  ggtitle("Pi per population - pink") +
  scale_fill_manual(values=c("red3", "orange", "cornflowerblue")) +
  # ggtitle("Pi per population - only new GVs phase III - pink") +
  # ylim(c(0.015,0.045)) +
  # ylim(c(0,0.02)) + # this one is for new GVs
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) + 
  # ggsave(paste0("plots/22_3_pi_pertreatment_total_pink.png"), width = 8, height = 7, dpi = 400)
  ggsave(paste0("plots/22_3_pi_pertreatment_newGV_pink.png"), width = 8, height = 7, dpi = 400)


test_table_pi<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>%
  filter(time %in% c(0, 53)) %>% 
  # filter(VarFreq>15) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  mutate(pi=(VarFreq*(100-VarFreq))/(((100-1)*100)/2)) %>% 
  # filter(time_groupGV=="phaseIII") %>%
  group_by(treatment, time, pop_ID) %>% 
  dplyr::summarize(N_GV=n(), total_pi=mean(pi)) %>% 
  ungroup() %>% 
  mutate(treatment=factor(treatment, 
                          levels=c("Allo_B", "Allo_T", "Para_B", "Para_T", "LM", "SYM", "P3P"))) %>%
  filter(treatment!="P3P") %>% 
  data.frame()

head(test_table_pi)
model<-lm(total_pi~treatment, test_table_pi) # to be run with phylogenetic control
summary(model)
anova(model)



####
# populations structure - Dxy

# Dxy per treatment:

# red populations:
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))


dxy_table_red<-data.frame(N_GV=c(), total_dxy=c(), pop1=c(), pop2=c(), treatment=c())
treatment_comparisons<-data.frame(
  treatment1=c("Allo_T", "Allo_T", "Allo_B", "Para_T", "Para_T", "Para_B", "LM", "LM", "LM", "SYM", "SYM", "SYM"), 
  treatment2=c("Allo_T", "Allo_B", "Allo_B", "Para_T", "Para_B", "Para_B", "LM", "Allo_T", "Allo_B", "SYM", "Allo_T", "Allo_B"))
subtable_phaseIII_table_total_red<-phaseIII_table_total_red %>%
  filter(time_groupGV!="fix_parental") %>% 
  #filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(53)) %>% 
  #filter(VarFreq>15) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:69) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  select(-treatment, -time)
# pairs parapatric populations:
population_para<-as.vector(unique(phaseIII_table_total_red$population[phaseIII_table_total_red$treatment=="Para_T"]))
treatment_used<-"Para_Pairs"
for (i in c(1:length(population_para))){
  pop1<-paste0(population_para[i], "_", "Para_T")
  pop2<-paste0(population_para[i], "_", "Para_B")
  dxy_table_red<-subtable_phaseIII_table_total_red %>% 
    filter(pop_ID %in% c(pop1, pop2)) %>% 
    spread(pop_ID, VarFreq, fill=0) %>% 
    mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
    mutate(dxy=((1-q_pop1)*q_pop2)+((1-q_pop2)*q_pop1)) %>% 
    ungroup() %>% 
    dplyr::summarize(N_GV=n(), total_dxy=mean(dxy)) %>% 
    mutate(pop1, pop2, treatment=treatment_used) %>% 
    rbind(dxy_table_red)
}
# Other comparisons:
for(line in seq(1,dim(treatment_comparisons)[1])){
  treatment1<-as.vector(treatment_comparisons[line,1])
  treatment2<-as.vector(treatment_comparisons[line,2])
  print(treatment1)
  print(treatment2)
  treatment_used<-paste0(treatment1, "__", treatment2)
  group_1<-unique(phaseIII_table_total_red$pop_treatment[phaseIII_table_total_red$treatment==treatment1])
  group_2<-unique(phaseIII_table_total_red$pop_treatment[phaseIII_table_total_red$treatment==treatment2])
  for (pop1_index in c(1:length(group_1))){
    for (pop2_index in c(1:length(group_2))){
      if (treatment1==treatment2){
        if(pop2_index>pop1_index){
          print(pop1_index)
          print(pop2_index)
          pop1<-group_1[pop1_index]
          pop2<-group_2[pop2_index]
          dxy_table_red<-subtable_phaseIII_table_total_red %>%
            filter(pop_ID %in% c(pop1, pop2)) %>% 
            spread(pop_ID, VarFreq, fill=0) %>% 
            mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
            mutate(dxy=((1-q_pop1)*q_pop2)+((1-q_pop2)*q_pop1)) %>% 
            ungroup() %>% 
            dplyr::summarize(N_GV=n(), total_dxy=mean(dxy)) %>% 
            mutate(pop1, pop2, treatment=treatment_used) %>% 
            rbind(dxy_table_red)
        }
      } else{
        print(pop1_index)
        print(pop2_index)
        pop1<-group_1[pop1_index]
        pop2<-group_2[pop2_index]
        dxy_table_red<-subtable_phaseIII_table_total_red %>%
          filter(pop_ID %in% c(pop1, pop2)) %>% 
          spread(pop_ID, VarFreq, fill=0) %>% 
          mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
          mutate(dxy=((1-q_pop1)*q_pop2)+((1-q_pop2)*q_pop1)) %>% 
          ungroup() %>% 
          dplyr::summarize(N_GV=n(), total_dxy=mean(dxy)) %>% 
          mutate(pop1, pop2, treatment=treatment_used) %>% 
          rbind(dxy_table_red)
      }
    }
  }
}

dxy_table_red$treatment<-factor(dxy_table_red$treatment, levels=c("Allo_B__Allo_B", "Allo_T__Allo_T", "Allo_T__Allo_B", "Para_T__Para_T", "Para_B__Para_B", "Para_T__Para_B", "Para_Pairs", "LM__LM", "LM__Allo_B", "LM__Allo_T", "SYM__SYM", "SYM__Allo_B", "SYM__Allo_T"))

write.table(dxy_table_red, "dxy_table_all_red.txt", quote = F, row.names=F)
dxy_table_red<-read.table("dxy_table_red.txt", T)


dxy_table_red %>%
  ggplot(aes(treatment, total_dxy)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.2) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        #axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+ 
  # ggsave(paste0("plots/22_3_dxy_total_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
  ggsave(paste0("plots/22_3_dxy_ancestralGV_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_dxy_allo_total_fixedVG_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_dxy_allo_total_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_dxy_allo_ancestralGV_phaseIII_red.png"), width = 8, height = 7, dpi = 400)

head(dxy_table_red)

library(lme4)
model<-lm(total_dxy~treatment, data=data.frame(dxy_table_red))
summary(model)  
anova(model)

# pink populations:
phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))


dxy_table_pink<-data.frame(N_GV=c(), total_dxy=c(), pop1=c(), pop2=c(), treatment=c())
treatment_comparisons<-data.frame(
  treatment1=c("Allo_T", "Allo_T", "Allo_B", "Para_T", "Para_T", "Para_B", "LM", "LM", "LM", "SYM", "SYM", "SYM"), 
  treatment2=c("Allo_T", "Allo_B", "Allo_B", "Para_T", "Para_B", "Para_B", "LM", "Allo_T", "Allo_B", "SYM", "Allo_T", "Allo_B"))
subtable_phaseIII_table_total_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(53)) %>% 
  #filter(VarFreq>15) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:69) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  select(-treatment, -time)
# pairs parapatric populations:
population_para<-as.vector(unique(phaseIII_table_total_pink$population[phaseIII_table_total_pink$treatment=="Para_T"]))
treatment_used<-"Para_Pairs"
for (i in c(1:length(population_para))){
  pop1<-paste0(population_para[i], "_", "Para_T")
  pop2<-paste0(population_para[i], "_", "Para_B")
  dxy_table_pink<-subtable_phaseIII_table_total_pink %>% 
    filter(pop_ID %in% c(pop1, pop2)) %>% 
    spread(pop_ID, VarFreq, fill=0) %>% 
    mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
    mutate(dxy=((1-q_pop1)*q_pop2)+((1-q_pop2)*q_pop1)) %>% 
    ungroup() %>% 
    dplyr::summarize(N_GV=n(), total_dxy=mean(dxy)) %>% 
    mutate(pop1, pop2, treatment=treatment_used) %>% 
    rbind(dxy_table_pink)
}
# Other comparisons:
for(line in seq(1,dim(treatment_comparisons)[1])){
  treatment1<-as.vector(treatment_comparisons[line,1])
  treatment2<-as.vector(treatment_comparisons[line,2])
  print(treatment1)
  print(treatment2)
  treatment_used<-paste0(treatment1, "__", treatment2)
  group_1<-unique(phaseIII_table_total_pink$pop_treatment[phaseIII_table_total_pink$treatment==treatment1])
  group_2<-unique(phaseIII_table_total_pink$pop_treatment[phaseIII_table_total_pink$treatment==treatment2])
  for (pop1_index in c(1:length(group_1))){
    for (pop2_index in c(1:length(group_2))){
      if (treatment1==treatment2){
        if(pop2_index>pop1_index){
          print(pop1_index)
          print(pop2_index)
          pop1<-group_1[pop1_index]
          pop2<-group_2[pop2_index]
          dxy_table_pink<-subtable_phaseIII_table_total_pink %>%
            filter(pop_ID %in% c(pop1, pop2)) %>% 
            spread(pop_ID, VarFreq, fill=0) %>% 
            mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
            mutate(dxy=((1-q_pop1)*q_pop2)+((1-q_pop2)*q_pop1)) %>% 
            ungroup() %>% 
            dplyr::summarize(N_GV=n(), total_dxy=mean(dxy)) %>% 
            mutate(pop1, pop2, treatment=treatment_used) %>% 
            rbind(dxy_table_pink)
        }
      } else{
        print(pop1_index)
        print(pop2_index)
        pop1<-group_1[pop1_index]
        pop2<-group_2[pop2_index]
        dxy_table_pink<-subtable_phaseIII_table_total_pink %>%
          filter(pop_ID %in% c(pop1, pop2)) %>% 
          spread(pop_ID, VarFreq, fill=0) %>% 
          mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
          mutate(dxy=((1-q_pop1)*q_pop2)+((1-q_pop2)*q_pop1)) %>% 
          ungroup() %>% 
          dplyr::summarize(N_GV=n(), total_dxy=mean(dxy)) %>% 
          mutate(pop1, pop2, treatment=treatment_used) %>% 
          rbind(dxy_table_pink)
      }
    }
  }
}

dxy_table_pink$treatment<-factor(dxy_table_pink$treatment, levels=c("Allo_B__Allo_B", "Allo_T__Allo_T", "Allo_T__Allo_B", "Para_T__Para_T", "Para_B__Para_B", "Para_T__Para_B", "Para_Pairs", "SYM__SYM", "SYM__Allo_B", "SYM__Allo_T", "LM__LM", "LM__Allo_B", "LM__Allo_T"))

dxy_table_pink %>% 
  ggplot(aes(treatment, total_dxy)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.2) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        #axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) + 
  # ggsave(paste0("plots/22_3_dxy_total_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
  ggsave(paste0("plots/22_3_dxy_ancestralGV_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_dxy_allo_total_fixedVG_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_dxy_allo_total_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_dxy_allo_ancestralGV_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)


### Population Structure - Fst
# red populations:
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

fst_table_red<-data.frame(N_GV=c(), total_fst=c(), pop1=c(), pop2=c(), treatment=c())
treatment_comparisons<-data.frame(
  treatment1=c("Allo_T", "Allo_T", "Allo_B", "Para_T", "Para_T", "Para_B", "LM", "LM", "LM", "SYM", "SYM", "SYM"), 
  treatment2=c("Allo_T", "Allo_B", "Allo_B", "Para_T", "Para_B", "Para_B", "LM", "Allo_T", "Allo_B", "SYM", "Allo_T", "Allo_B"))

subtable_phaseIII_table_total_red<-
  phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(53)) %>% 
  # filter(VarFreq>15) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:69) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  select(-treatment, -time)

population_para<-as.vector(unique(phaseIII_table_total_red$population[phaseIII_table_total_red$treatment=="Para_T"]))
treatment_used<-"Para_Pairs"
for (i in c(1:length(population_para))){
  pop1<-paste0(population_para[i], "_", "Para_T")
  pop2<-paste0(population_para[i], "_", "Para_B")
  fst_table_red<-subtable_phaseIII_table_total_red %>%
    filter(pop_ID %in% c(pop1, pop2)) %>% 
    spread(pop_ID, VarFreq, fill=0) %>% 
    mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100, 
           pi_pop1=((1-q_pop1)*q_pop1)/(((100-1)*100)/2), 
           pi_pop2=((1-q_pop2)*q_pop2)/(((100-1)*100)/2), 
           mean_q=mean(c(q_pop1, q_pop2)), 
           pi_total=((1-mean_q)*mean_q)/(((100-1)*100)/2), 
           fst=(pi_total-(mean(c(pi_pop1,pi_pop2))))/pi_total, 
           fst=ifelse(is.nan(fst), 0, fst)) %>% 
    ungroup() %>% 
    dplyr::summarize(N_GV=n(), total_fst=mean(fst)) %>% 
    mutate(pop1, pop2, treatment=treatment_used) %>% 
    rbind(fst_table_red)
}
# Other comparisons:
for(line in seq(1,dim(treatment_comparisons)[1])){
  treatment1<-as.vector(treatment_comparisons[line,1])
  treatment2<-as.vector(treatment_comparisons[line,2])
  print(treatment1)
  print(treatment2)
  treatment_used<-paste0(treatment1, "__", treatment2)
  group_1<-unique(phaseIII_table_total_red$pop_treatment[phaseIII_table_total_red$treatment==treatment1])
  group_2<-unique(phaseIII_table_total_red$pop_treatment[phaseIII_table_total_red$treatment==treatment2])
  for (pop1_index in c(1:length(group_1))){
    for (pop2_index in c(1:length(group_2))){
      if (treatment1==treatment2){
        if(pop2_index>pop1_index){
          print(pop1_index)
          print(pop2_index)
          pop1<-group_1[pop1_index]
          pop2<-group_2[pop2_index]
          fst_table_red<-subtable_phaseIII_table_total_red %>%
            filter(pop_ID %in% c(pop1, pop2)) %>% 
            spread(pop_ID, VarFreq, fill=0) %>% 
            mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
            mutate(pi_pop1=((1-q_pop1)*q_pop1)/(((100-1)*100)/2), 
                   pi_pop2=((1-q_pop2)*q_pop2)/(((100-1)*100)/2), 
                   mean_q=mean(c(q_pop1, q_pop2)), 
                   pi_total=((1-mean_q)*mean_q)/(((100-1)*100)/2), 
                   fst=(pi_total-(mean(c(pi_pop1,pi_pop2))))/pi_total) %>% 
            mutate(fst=ifelse(is.nan(fst), 0, fst)) %>% 
            ungroup() %>% 
            dplyr::summarize(N_GV=n(), total_fst=mean(fst)) %>% 
            mutate(pop1, pop2, treatment=treatment_used) %>% 
            rbind(fst_table_red)
        }
      } else{
        print(pop1_index)
        print(pop2_index)
        pop1<-group_1[pop1_index]
        pop2<-group_2[pop2_index]
        fst_table_red<-subtable_phaseIII_table_total_red %>%
          filter(pop_ID %in% c(pop1, pop2)) %>% 
          spread(pop_ID, VarFreq, fill=0) %>% 
          mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
          mutate(pi_pop1=((1-q_pop1)*q_pop1)/(((100-1)*100)/2), 
                 pi_pop2=((1-q_pop2)*q_pop2)/(((100-1)*100)/2), 
                 mean_q=mean(c(q_pop1, q_pop2)), 
                 pi_total=((1-mean_q)*mean_q)/(((100-1)*100)/2), 
                 fst=(pi_total-(mean(c(pi_pop1,pi_pop2))))/pi_total) %>% 
          mutate(fst=ifelse(is.nan(fst), 0, fst)) %>% 
          ungroup() %>% 
          dplyr::summarize(N_GV=n(), total_fst=mean(fst)) %>% 
          mutate(pop1, pop2, treatment=treatment_used) %>% 
          rbind(fst_table_red)
      }
    }
  }
}


fst_table_red$treatment<-factor(fst_table_red$treatment, levels=c("Allo_B__Allo_B", "Allo_T__Allo_T", "Allo_T__Allo_B", "Para_T__Para_T", "Para_B__Para_B", "Para_T__Para_B", "Para_Pairs", "SYM__SYM", "SYM__Allo_B", "SYM__Allo_T", "LM__LM", "LM__Allo_B", "LM__Allo_T"))

fst_table_red %>% 
  ggplot(aes(treatment, total_fst)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.2) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        #axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) + 
  ggsave(paste0("plots/22_3_fst_total_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
  # ggsave(paste0("plots/22_3_fst_ancestralGV_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_fst_allo_total_fixedVG_phaseIII_red.png"),width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_fst_allo_total_phaseIII_red.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_fst_allo_ancestralGV_phaseIII_red.png"), width = 8, height = 7, dpi = 400)




# pink populations:

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

fst_table_pink<-data.frame(N_GV=c(), total_fst=c(), pop1=c(), pop2=c(), treatment=c())
treatment_comparisons<-data.frame(
  treatment1=c("Allo_T", "Allo_T", "Allo_B", "Para_T", "Para_T", "Para_B", "LM", "LM", "LM", "SYM", "SYM", "SYM"), 
  treatment2=c("Allo_T", "Allo_B", "Allo_B", "Para_T", "Para_B", "Para_B", "LM", "Allo_T", "Allo_B", "SYM", "Allo_T", "Allo_B"))

subtable_phaseIII_table_total_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(53)) %>% 
  # filter(VarFreq>15) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>%  
  gather("time_pop_ID", "VarFreq", 5:69) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  select(-treatment, -time)

population_para<-as.vector(unique(phaseIII_table_total_pink$population[phaseIII_table_total_pink$treatment=="Para_T"]))
treatment_used<-"Para_Pairs"
for (i in c(1:length(population_para))){
  pop1<-paste0(population_para[i], "_", "Para_T")
  pop2<-paste0(population_para[i], "_", "Para_B")
  fst_table_pink<-subtable_phaseIII_table_total_pink %>%
    filter(pop_ID %in% c(pop1, pop2)) %>% 
    spread(pop_ID, VarFreq, fill=0) %>% 
    mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100, 
           pi_pop1=((1-q_pop1)*q_pop1)/(((100-1)*100)/2), 
           pi_pop2=((1-q_pop2)*q_pop2)/(((100-1)*100)/2), 
           mean_q=mean(c(q_pop1, q_pop2)), 
           pi_total=((1-mean_q)*mean_q)/(((100-1)*100)/2), 
           fst=(pi_total-(mean(c(pi_pop1,pi_pop2))))/pi_total, 
           fst=ifelse(is.nan(fst), 0, fst)) %>% 
    ungroup() %>% 
    dplyr::summarize(N_GV=n(), total_fst=mean(fst)) %>% 
    mutate(pop1, pop2, treatment=treatment_used) %>% 
    rbind(fst_table_pink)
}
# Other comparisons:
for(line in seq(1,dim(treatment_comparisons)[1])){
  treatment1<-as.vector(treatment_comparisons[line,1])
  treatment2<-as.vector(treatment_comparisons[line,2])
  print(treatment1)
  print(treatment2)
  treatment_used<-paste0(treatment1, "__", treatment2)
  group_1<-unique(phaseIII_table_total_pink$pop_treatment[phaseIII_table_total_pink$treatment==treatment1])
  group_2<-unique(phaseIII_table_total_pink$pop_treatment[phaseIII_table_total_pink$treatment==treatment2])
  for (pop1_index in c(1:length(group_1))){
    for (pop2_index in c(1:length(group_2))){
      if (treatment1==treatment2){
        if(pop2_index>pop1_index){
          print(pop1_index)
          print(pop2_index)
          pop1<-group_1[pop1_index]
          pop2<-group_2[pop2_index]
          fst_table_pink<-subtable_phaseIII_table_total_pink %>%
            filter(pop_ID %in% c(pop1, pop2)) %>% 
            spread(pop_ID, VarFreq, fill=0) %>% 
            mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
            mutate(pi_pop1=((1-q_pop1)*q_pop1)/(((100-1)*100)/2), 
                   pi_pop2=((1-q_pop2)*q_pop2)/(((100-1)*100)/2), 
                   mean_q=mean(c(q_pop1, q_pop2)), 
                   pi_total=((1-mean_q)*mean_q)/(((100-1)*100)/2), 
                   fst=(pi_total-(mean(c(pi_pop1,pi_pop2))))/pi_total) %>% 
            mutate(fst=ifelse(is.nan(fst), 0, fst)) %>% 
            ungroup() %>% 
            dplyr::summarize(N_GV=n(), total_fst=mean(fst)) %>% 
            mutate(pop1, pop2, treatment=treatment_used) %>% 
            rbind(fst_table_pink)
        }
      } else{
        print(pop1_index)
        print(pop2_index)
        pop1<-group_1[pop1_index]
        pop2<-group_2[pop2_index]
        fst_table_pink<-subtable_phaseIII_table_total_pink %>%
          filter(pop_ID %in% c(pop1, pop2)) %>% 
          spread(pop_ID, VarFreq, fill=0) %>% 
          mutate(q_pop1=get(pop1)/100, q_pop2=get(pop2)/100) %>% 
          mutate(pi_pop1=((1-q_pop1)*q_pop1)/(((100-1)*100)/2), 
                 pi_pop2=((1-q_pop2)*q_pop2)/(((100-1)*100)/2), 
                 mean_q=mean(c(q_pop1, q_pop2)), 
                 pi_total=((1-mean_q)*mean_q)/(((100-1)*100)/2), 
                 fst=(pi_total-(mean(c(pi_pop1,pi_pop2))))/pi_total) %>% 
          mutate(fst=ifelse(is.nan(fst), 0, fst)) %>% 
          ungroup() %>% 
          dplyr::summarize(N_GV=n(), total_fst=mean(fst)) %>% 
          mutate(pop1, pop2, treatment=treatment_used) %>% 
          rbind(fst_table_pink)
      }
    }
  }
}


fst_table_pink$treatment<-factor(fst_table_pink$treatment, levels=c("Allo_B__Allo_B", "Allo_T__Allo_T", "Allo_T__Allo_B", "Para_T__Para_T", "Para_B__Para_B", "Para_T__Para_B", "Para_Pairs", "SYM__SYM", "SYM__Allo_B", "SYM__Allo_T", "LM__LM", "LM__Allo_B", "LM__Allo_T"))

fst_table_pink %>% 
  ggplot(aes(treatment, total_fst)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.2) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        #axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave(paste0("plots/22_3_fst_total_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
  # ggsave(paste0("plots/22_3_fst_ancestralGV_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_fst_allo_total_fixedVG_phaseIII_pink.png"),width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_fst_allo_total_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)
# ggsave(paste0("plots/22_3_fst_allo_ancestralGV_phaseIII_pink.png"), width = 8, height = 7, dpi = 400)





#####
# Number of populations with especific genetic variants:
# red:
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV) %>% 
  filter(VarFreq>0) %>% 
  dplyr::summarize(N_population=n(), VarFreq=mean(VarFreq)) %>% 
  #filter(N_population>1 & VarFreq>80 & time_groupGV=="phaseIII")
  ggplot(aes(N_population, VarFreq, colour=time_groupGV)) + 
  geom_point(size=2, alpha=0.5) +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(time_groupGV ~ .) #+ 
  # ggsave(paste0("plots/22_4_Number_pop_perGV_red.png"), width = 8, height = 10, dpi = 400)

b<-seq(0,1,0.05)
cuts <- c(-Inf, b[-1]-diff(b)/2, Inf)

phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  mutate(VarFreq=VarFreq/100) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>%
  filter(time_groupGV=="phaseIII") %>% 
  # # table:
  # filter(as.vector(VarFreq) > 0) %>%
  # group_by(time_pop_ID) %>%
  # dplyr::summarize(N_GV=n(),
  #                  VarFreq_25=sum((VarFreq>0.25)*1),
  #                  proportion_var25=VarFreq_25/N_GV)  %>%
  # summary()
  # # plot:
  group_by(Chrom, Position, GV_code, time_groupGV) %>% 
  mutate(VarFreq_cat=cut(VarFreq, breaks=cuts, labels=b)) %>% 
  group_by(time_pop_ID, VarFreq_cat) %>%
  dplyr::summarize(N_GV=n()) %>% 
  filter(as.vector(VarFreq_cat) > 0) %>%
  separate(time_pop_ID, c("treatment", "time", "pop"), sep="__") %>%
  mutate(treatment=factor(treatment, levels=c("Allo_T", "Allo_B", "Para_T", 
                                    "Para_B", "LM", "SYM"))) %>%
  ggplot(aes(VarFreq_cat, N_GV)) +
  geom_boxplot() +
  xlab("Allele Fq.") +
  ylab("Num. Genetic Variants") +
  ylab("Num. New Genetic Variants") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title = element_text(size=14, colour="black"), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+ 
  ggsave(paste0("plots/22_4_SFS_newGV_red.png"), width = 6, height = 3, dpi = 400)
  ggsave(paste0("plots/22_4_SFS_red.png"), width = 6, height = 3, dpi = 400)

phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  mutate(VarFreq=VarFreq/100) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>%
  filter(time_groupGV=="phaseIII") %>% 
  group_by(Chrom, Position, GV_code, time_groupGV) %>% 
  mutate(VarFreq_cat=cut(VarFreq, breaks=cuts, labels=b)) %>% 
  # filter(as.vector(VarFreq_cat) > 0) %>%
  group_by(time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), 
                   low_25=(sum((VarFreq<0.26)*1)/N_GV), 
                   low_30=(sum((VarFreq<0.3)*1)/N_GV), 
                   low_40=(sum((VarFreq<0.4)*1)/N_GV), 
                   high_85=(sum((VarFreq>0.85)*1)/N_GV), 
                   high_90=(sum((VarFreq>0.9)*1)/N_GV), 
                   Total_low25_plu85=low_25+high_85) %>% 
  # data:
  ungroup() %>%
  dplyr::summarize(N_GV=n(),
                   # min_25=min(low_25),
                   # max_25=max(low_25),
                   # min_85=min(high_85),
                   # max_85=max(high_85),
                   # min_total=min(Total_low25_plu85),
                   # max_total=max(Total_low25_plu85),
                   mean_25=mean(low_25),
                   sd_25=sd(low_25),
                   mean_85=mean(high_85),
                   sd_85=sd(high_85),
                   mean_total=mean(Total_low25_plu85),
                   sd_total=sd(Total_low25_plu85))
  # # plot:
  # ggplot(aes(low_25)) +
  # geom_histogram() +
  # geom_histogram(aes(high_85), fill="red") +
  # geom_histogram(aes(Total_low25_plu85), fill="blue", alpha=0.7) +
  # xlab("Proportion") +
  # ylab("Count") +
  # # xlim(c(0,1)) + 
  # theme_classic() + 
  # theme(panel.border = element_blank(), 
  #       panel.grid.major = element_line(colour = "grey90"), 
  #       panel.grid.minor = element_line(colour = "grey90"), 
  #       strip.background = element_blank(), 
  #       strip.text = element_blank(), 
  #       axis.title = element_text(size=14, colour="black"), 
  #       axis.line = element_line(colour = "black"), 
  #       axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
  #       axis.text.y = element_text(size=12, colour="black")) #+ 
  # ggsave(paste0("plots/22_4_SFS_NewGV_proportions_25_85_red.png"), width = 3, height = 3, dpi = 400)
  # ggsave(paste0("plots/22_4_SFS_proportions_25_85_red.png"), width = 3, height = 3, dpi = 400)

  
# SFS for fractions:
  
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53, 55, 57) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))




library(ggrepel)
library(wesanderson)

treatment<-c()
populationID<-c()
Pool<-c()
Top<-c()
Bottom<-c()


treatment_used="LM"

populations<-phaseIII_table_total_red %>% 
  filter(treatment==treatment_used) %>% 
  ungroup() %>% 
  select(population) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()


for (population_used in populations){
  print(population_used)
  # if(!(population_used %in% c("C10"))) { # for SYM populations
  temporal<-phaseIII_table_total_red %>% 
    mutate(VarFreq=VarFreq/100) %>% 
    filter(treatment==treatment_used) %>% 
    filter(time_groupGV!="fix_parental") %>% 
    filter(population==population_used) %>% 
    mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
    ungroup() %>% 
    select(time, GV_code, VarFreq) %>% 
    spread(time, VarFreq, fill=0) %>% 
    # filter(`53`!=0) %>% 
    mutate(namesID=if_else(`53`>0.20, as.character(round(`53`, digits = 2)), "")) %>% 
    mutate(Pool=`53`, 
           Top=`55`, 
           Bottom=`57`) 
  treatment<-c(treatment, rep(treatment_used, dim(temporal)[1]))
  populationID<-c(populationID, rep(population_used, dim(temporal)[1]))
  Pool<-c(Pool, temporal$Pool)
  Top<-c(Top, temporal$Top)
  Bottom<-c(Bottom, temporal$Bottom)
  # temporal %>%
  #   ggplot(aes(Bottom, Top, colour=Pool)) +
  #   geom_point(size=2)+
  #   geom_abline(slope=1, intercept=0, colour='grey80') +
  #   geom_text_repel(aes(label = namesID),
  #                   box.padding   = 0.35,
  #                   point.padding = 0.5,
  #                   segment.color = 'grey40',
  #                   colour="black") +
  #   ggtitle(paste0("Red ", treatment_used, " Pop=", population_used))+
  #   scale_colour_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous")) +
  #   scale_x_continuous(breaks=seq(0,1,0.2))+
  #   scale_y_continuous(breaks=seq(0,1,0.2))+
  #   theme_classic() +
  #   theme(panel.border = element_blank(),
  #         panel.grid.major = element_line(colour = "grey90"),
  #         panel.grid.minor = element_line(colour = "grey90"),
  #         axis.line = element_line(colour = "black"),
  #         legend.position="right",
  #         axis.text.x = element_text(size=12, colour="black"),
  #         axis.text.y = element_text(size=12, colour="black"))+
  #   ggsave(paste0("plots/26_comparison_fractions_", treatment_used, "_pop", population_used, "_red.png"), width = 7, height = 6, dpi = 400)
  # }
}


data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>%
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
  #           N_GV_Top=sum(diff_top_bottom>0.15), 
  #           N_GV_Bottom=sum(diff_top_bottom<(-0.15))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
            N_GV_Top=sum(diff_top_bottom>0.15), 
            N_GV_Bottom=sum(diff_top_bottom<(-0.15))*(-1),
            prop_top=N_GV_Top/N_GV, 
            prop_bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(value_top=paste(N_GV_Top, prop_top, sep="__"), 
         value_bottom=paste(N_GV_Bottom, prop_bottom, sep="__")) %>% 
  select(treatment, populationID, value_top, value_bottom) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  # filter(group!="N_GV") %>% 
  # filter(group!="N_GV_NC") %>% 
  # filter(group %in% c("prop_top", "prop_bottom")) %>%
  # filter(group %in% c("N_GV_Top", "N_GV_Bottom")) %>%
  ggplot(aes(populationID, Porportion, fill=group))+
  # geom_bar(stat="identity",position=position_dodge())+
  geom_bar(stat="identity")+
  geom_text(aes(label = abs(N_GV), 
                vjust = ifelse(Porportion >= 0, (-0.5), 1))) +
  scale_fill_manual(values=c("red3", "cornflowerblue"))+
  xlab("Population") +
  # ggplot(aes(group, NGenVar, fill=group))+
  # geom_boxplot()+
  # geom_jitter(alpha=0.9) +
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=90,size=12,colour="black",hjust= 1,vjust=0.5),
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave("plots/27_dominantFraction_red.png", width = 7, height = 5, dpi = 400)

data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>%
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
  #           N_GV_Top=sum(diff_top_bottom>0.15), 
  #           N_GV_Bottom=sum(diff_top_bottom<(-0.15))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
            N_GV_Top=sum(diff_top_bottom>0.15), 
            N_GV_Bottom=sum(diff_top_bottom<(-0.15)),
            prop_top=N_GV_Top/N_GV, 
            prop_bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Top=paste(N_GV_Top, prop_top, sep="__"), 
         diff_Bot=paste(N_GV_Bottom, prop_bottom, sep="__")) %>% 
  select(treatment, populationID, diff_Top, diff_Bot) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  ggplot(aes(group, Porportion, fill=group))+
  geom_boxplot()+
  scale_fill_manual(values=c("red3", "cornflowerblue"))+
  facet_grid(. ~ treatment, space="free", scale="free")+
  ylim(value=c(0,0.35)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("plots/27_dominantFraction_boxplot_red.png", width = 3, height = 4, dpi = 400)




data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_pool_top>(-0.15) & diff_pool_top<0.15)),
  #           N_GV_Pool=sum(diff_pool_top>0.15), 
  #           N_GV_Top=sum(diff_pool_top<(-0.15))*(-1),
  #           prop_Pool=N_GV_Pool/N_GV, 
  #           prop_Top=((N_GV_Top/N_GV))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_top>(-0.20) & diff_pool_top<0.20)),
            N_GV_Pool=sum(diff_pool_top>0.20), 
            N_GV_Top=sum(diff_pool_top<(-0.20))*(-1),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Top=((N_GV_Top/N_GV))) %>% 
  mutate(value_Pool=paste(N_GV_Pool, prop_Pool, sep="__"), 
         value_Top=paste(N_GV_Top, prop_Top, sep="__")) %>% 
  select(treatment, populationID, value_Pool, value_Top) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  # filter(group!="N_GV") %>% 
  # filter(group!="N_GV_NC") %>% 
  # filter(group %in% c("prop_top", "prop_bottom")) %>%
  # filter(group %in% c("N_GV_Pool", "N_GV_Top")) %>%
  ggplot(aes(populationID, Porportion, fill=group))+
  # geom_bar(stat="identity",position=position_dodge())+
  geom_bar(stat="identity")+
  geom_text(aes(label = abs(N_GV), 
                vjust = ifelse(Porportion >= 0, (-0.5), 1))) +
  scale_fill_manual(values=c("orange", "cornflowerblue"))+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  xlab("Population") +
  # ggplot(aes(group, NGenVar, fill=group))+
  # geom_boxplot()+
  # geom_jitter(alpha=0.9) +
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=90,size=12,colour="black",hjust= 1,vjust=0.5),
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("plots/27_dominantFraction_pool_top_red20.png", width = 7, height = 5, dpi = 400)


data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_pool_bottom>(-0.15) & diff_pool_bottom<0.15)),
  #           N_GV_Pool=sum(diff_pool_bottom>0.15), 
  #           N_GV_Bottom=sum(diff_pool_bottom<(-0.15))*(-1),
  #           prop_Pool=N_GV_Pool/N_GV, 
  #           prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_bottom>(-0.20) & diff_pool_bottom<0.20)),
            N_GV_Pool=sum(diff_pool_bottom>0.20), 
            N_GV_Bottom=sum(diff_pool_bottom<(-0.20))*(-1),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(value_Pool=paste(N_GV_Pool, prop_Pool, sep="__"), 
         value_Bottom=paste(N_GV_Bottom, prop_Bottom, sep="__")) %>% 
  select(treatment, populationID, value_Pool, value_Bottom) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  # filter(group!="N_GV") %>% 
  # filter(group!="N_GV_NC") %>% 
  # filter(group %in% c("prop_top", "prop_bottom")) %>%
  # filter(group %in% c("N_GV_Pool", "N_GV_Bottom")) %>%
  ggplot(aes(populationID, Porportion, fill=group))+
  # geom_bar(stat="identity",position=position_dodge())+
  geom_bar(stat="identity")+
  geom_text(aes(label = abs(N_GV), 
                vjust = ifelse(Porportion >= 0, (-0.5), 1))) +
  scale_fill_manual(values=c("red3", "orange"))+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  xlab("Population") +
  # ggplot(aes(group, NGenVar, fill=group))+
  # geom_boxplot()+
  # geom_jitter(alpha=0.9) +
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=90,size=12,colour="black",hjust= 1,vjust=0.5),
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("plots/27_dominantFraction_pool_bottom_red20.png", width = 7, height = 5, dpi = 400)



data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_top>(-0.20) & diff_pool_top<0.20)),
            N_GV_Pool=sum(diff_pool_top>0.20), 
            N_GV_Top=sum(diff_pool_top<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Top=((N_GV_Top/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Top=paste(N_GV_Top, prop_Top, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Top) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  ggplot(aes(group, Porportion, fill=group))+
  geom_boxplot()+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  scale_fill_manual(values=c("orange", "cornflowerblue"))+
  facet_grid(. ~ treatment, space="free", scale="free")+
  ylim(value=c(0,0.35)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black")) #+
# ggsave("plots/27_dominantFraction_boxplot_pool_top_red.png", width = 3, height = 4, dpi = 400)


data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_bottom>(-0.20) & diff_pool_bottom<0.20)),
            N_GV_Pool=sum(diff_pool_bottom>0.20), 
            N_GV_Bottom=sum(diff_pool_bottom<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Bot=paste(N_GV_Bottom, prop_Bottom, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Bot) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  mutate(group=factor(group, levels=c("diff_Pool", "diff_Bot"))) %>% 
  ggplot(aes(group, Porportion, fill=group))+
  geom_boxplot()+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  scale_fill_manual(values=c("orange", "red3"))+
  facet_grid(. ~ treatment, space="free", scale="free")+
  ylim(value=c(0,0.35)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("plots/27_dominantFraction_boxplot_pool_bottom_red.png", width = 3, height = 4, dpi = 400)


# Statistical test:

# red - top - bottom comparison
table_dominanFraction_top_bottom_red<-data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>%
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
            N_GV_Top=sum(diff_top_bottom>0.15), 
            N_GV_Bottom=sum(diff_top_bottom<(-0.15)),
            prop_top=N_GV_Top/N_GV, 
            prop_bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Top=paste(N_GV_Top, prop_top, sep="__"), 
         diff_Bot=paste(N_GV_Bottom, prop_bottom, sep="__")) %>% 
  select(treatment, populationID, diff_Top, diff_Bot) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  select(-populationID) %>%
  data.frame()

mod <- glm(Porportion~treatment/group,family=quasibinomial,table_dominanFraction_top_bottom_red)
summary(mod)

mod <- glm(Porportion~treatment + group + treatment*group,family=quasibinomial,table_dominanFraction_top_bottom_red)
summary(mod)

# red - Pool - top bottom comparison
table_dominanFraction_pool_top_red<-data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_top>(-0.20) & diff_pool_top<0.20)),
            N_GV_Pool=sum(diff_pool_top>0.20), 
            N_GV_Top=sum(diff_pool_top<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Top=((N_GV_Top/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Top=paste(N_GV_Top, prop_Top, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Top) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  data.frame()
 
mod <- glm(Porportion~treatment/group,family=quasibinomial,table_dominanFraction_pool_top_red)
summary(mod) 

# red - Pool - Bottom bottom comparison
table_dominanFraction_pool_bottom_red<-data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_bottom>(-0.20) & diff_pool_bottom<0.20)),
            N_GV_Pool=sum(diff_pool_bottom>0.20), 
            N_GV_Bottom=sum(diff_pool_bottom<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Bot=paste(N_GV_Bottom, prop_Bottom, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Bot) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  mutate(group=factor(group, levels=c("diff_Pool", "diff_Bot"))) %>% 
  data.frame()

mod <- glm(Porportion~treatment/group,family=quasibinomial,table_dominanFraction_pool_bottom_red)
summary(mod) 




# pink:

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV) %>% 
  filter(VarFreq>0) %>% 
  dplyr::summarize(N_population=n(), VarFreq=mean(VarFreq)) %>% 
  #filter(N_population>1 & VarFreq>80 & time_groupGV=="phaseIII")
  ggplot(aes(N_population, VarFreq, colour=time_groupGV)) + 
  geom_point(size=2, alpha=0.5) +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  facet_grid(time_groupGV ~ .) #+ 
  ggsave(paste0("plots/22_4_Number_pop_perGV_pink.png"), width = 8, height = 10, dpi = 400)

b<-seq(0,1,0.05)
cuts <- c(-Inf, b[-1]-diff(b)/2, Inf)

phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  mutate(VarFreq=VarFreq/100) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>%
  filter(time_groupGV=="phaseIII") %>% 
  # # table:
  # filter(as.vector(VarFreq) > 0) %>%
  # group_by(time_pop_ID) %>%
  # dplyr::summarize(N_GV=n(),
  #                  VarFreq_25=sum((VarFreq>0.25)*1),
  #                  proportion_var25=VarFreq_25/N_GV)  %>%
  # summary()
  # # plot:
  group_by(Chrom, Position, GV_code, time_groupGV) %>% 
  mutate(VarFreq_cat=cut(VarFreq, breaks=cuts, labels=b)) %>% 
  group_by(time_pop_ID, VarFreq_cat) %>%
  dplyr::summarize(N_GV=n()) %>% 
  filter(as.vector(VarFreq_cat) > 0) %>%
  separate(time_pop_ID, c("treatment", "time", "pop"), sep="__") %>%
  mutate(treatment=factor(treatment, levels=c("Allo_T", "Allo_B", "Para_T", 
                                              "Para_B", "LM", "SYM"))) %>%
  ggplot(aes(VarFreq_cat, N_GV)) +
  geom_boxplot() +
  xlab("Allele Fq.") +
  ylab("Num. New Genetic Variants") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title = element_text(size=14, colour="black"), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) + 
  ggsave(paste0("plots/22_4_SFS_newGV_pink.png"), width = 6, height = 3, dpi = 400)
  ggsave(paste0("plots/22_4_SFS_pink.png"), width = 6, height = 3, dpi = 400)

phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>%
  mutate(VarFreq=VarFreq/100) %>%
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>%
  group_by(Chrom, Position, GV_code, time_groupGV) %>% 
  mutate(VarFreq_cat=cut(VarFreq, breaks=cuts, labels=b)) %>% 
  group_by(time_pop_ID) %>% 
  dplyr::summarize(N_GV=n(), 
                   low_25=(sum((VarFreq<0.26)*1)/N_GV), 
                   low_30=(sum((VarFreq<0.3)*1)/N_GV), 
                   low_40=(sum((VarFreq<0.4)*1)/N_GV), 
                   high_85=(sum((VarFreq>0.85)*1)/N_GV), 
                   high_90=(sum((VarFreq>0.9)*1)/N_GV), 
                   Total_low25_plu85=low_25+high_85) %>% 
  # table:
  ungroup() %>%
  dplyr::summarize(N_GV=n(), 
                   # min_25=min(low_25), 
                   # max_25=max(low_25), 
                   # min_85=min(high_85), 
                   # max_85=max(high_85), 
                   # min_total=min(Total_low25_plu85), 
                   # max_total=max(Total_low25_plu85), 
                   mean_25=mean(low_25),
                   sd_25=sd(low_25), 
                   mean_85=mean(high_85),
                   sd_85=sd(high_85), 
                   mean_total=mean(Total_low25_plu85),
                   sd_total=sd(Total_low25_plu85))
  # plot:
  ggplot(aes(low_25)) +
  geom_histogram() +
  geom_histogram(aes(high_85), fill="red") +
  geom_histogram(aes(Total_low25_plu85), fill="blue", alpha=0.7) +
  xlab("Proportion") +
  ylab("Count") +
  xlim(c(0,1)) + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title = element_text(size=14, colour="black"), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) + 
  ggsave(paste0("plots/22_4_SFS_proportions_25_85_pink.png"), width = 3, height = 3, dpi = 400)




# this part is only to check individual GVs:
rbind(table_indels, table_snps) %>% 
  #select(GV_code, time, population, treatment, parental, VarFreq) %>% 
  filter(time %in% c((-15), (-16), (-10), (-5), 0, 53) & 
           parental=="P3R" &
           Position==3845360) -> temporal
#GV_code=="I_3845360_C") -> temporal


# SFS for fractions pink:

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53, 55, 57) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))



library(ggrepel)
library(wesanderson)

treatment<-c()
populationID<-c()
Pool<-c()
Top<-c()
Bottom<-c()

treatment_used="LM"

populations<-phaseIII_table_total_pink %>% 
  filter(treatment==treatment_used) %>% 
  ungroup() %>% 
  select(population) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()

for (population_used in populations){
  print(population_used)
  temporal<-phaseIII_table_total_pink %>% 
    mutate(VarFreq=VarFreq/100) %>% 
    filter(treatment==treatment_used) %>% 
    filter(time_groupGV!="fix_parental") %>% 
    filter(population==population_used) %>% 
    mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
    group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
    dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
    select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
    separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
    ungroup() %>% 
    select(time, GV_code, VarFreq) %>% 
    spread(time, VarFreq, fill=0) %>% 
    # filter(`53`!=0) %>% 
    mutate(namesID=if_else(`53`>0.20, as.character(round(`53`, digits = 2)), "")) %>% 
    mutate(Pool=`53`, 
           Top=`55`, 
           Bottom=`57`) 
  treatment<-c(treatment, rep(treatment_used, dim(temporal)[1]))
  populationID<-c(populationID, rep(population_used, dim(temporal)[1]))
  Pool<-c(Pool, temporal$Pool)
  Top<-c(Top, temporal$Top)
  Bottom<-c(Bottom, temporal$Bottom)
  # temporal %>%
  #   ggplot(aes(Bottom, Top, colour=Pool)) +
  #   geom_point(size=2)+
  #   geom_abline(slope=1, intercept=0, colour='grey80') +
  #   geom_text_repel(aes(label = namesID),
  #                   box.padding   = 0.35,
  #                   point.padding = 0.5,
  #                   segment.color = 'grey40',
  #                   colour="black") +
  #   ggtitle(paste0("pink ", treatment_used, " Pop=", population_used))+
  #   scale_colour_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous")) +
  #   scale_x_continuous(breaks=seq(0,1,0.2))+
  #   scale_y_continuous(breaks=seq(0,1,0.2))+
  #   theme_classic() +
  #   theme(panel.border = element_blank(),
  #         panel.grid.major = element_line(colour = "grey90"),
  #         panel.grid.minor = element_line(colour = "grey90"),
  #         axis.line = element_line(colour = "black"),
  #         legend.position="right",
  #         axis.text.x = element_text(size=12, colour="black"),
  #         axis.text.y = element_text(size=12, colour="black"))+
  #   ggsave(paste0("plots/26_comparison_fractions_", treatment_used, "_pop", population_used, "_pink.png"), width = 7, height = 6, dpi = 400)
  # }
}

data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>%
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
  #           N_GV_Top=sum(diff_top_bottom>0.15), 
  #           N_GV_Bottom=sum(diff_top_bottom<(-0.15))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
            N_GV_Top=sum(diff_top_bottom>0.15), 
            N_GV_Bottom=sum(diff_top_bottom<(-0.15))*(-1),
            prop_top=N_GV_Top/N_GV, 
            prop_bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(value_top=paste(N_GV_Top, prop_top, sep="__"), 
         value_bottom=paste(N_GV_Bottom, prop_bottom, sep="__")) %>% 
  select(treatment, populationID, value_top, value_bottom) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  # filter(group!="N_GV") %>% 
  # filter(group!="N_GV_NC") %>% 
  # filter(group %in% c("prop_top", "prop_bottom")) %>%
  # filter(group %in% c("N_GV_Top", "N_GV_Bottom")) %>%
  ggplot(aes(populationID, Porportion, fill=group))+
  # geom_bar(stat="identity",position=position_dodge())+
  geom_bar(stat="identity")+
  geom_text(aes(label = abs(N_GV), 
                vjust = ifelse(Porportion >= 0, (-0.5), 1))) +
  scale_fill_manual(values=c("red3", "cornflowerblue"))+
  xlab("Population") +
  # ggplot(aes(group, NGenVar, fill=group))+
  # geom_boxplot()+
  # geom_jitter(alpha=0.9) +
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=90,size=12,colour="black",hjust= 1,vjust=0.5),
        axis.text.y = element_text(size=12, colour="black"))#+
  ggsave("plots/27_dominantFraction_pink.png", width = 7, height = 5, dpi = 400)

data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>%
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
  #           N_GV_Top=sum(diff_top_bottom>0.15), 
  #           N_GV_Bottom=sum(diff_top_bottom<(-0.15))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
            N_GV_Top=sum(diff_top_bottom>0.15), 
            N_GV_Bottom=sum(diff_top_bottom<(-0.15)),
            prop_top=N_GV_Top/N_GV, 
            prop_bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Top=paste(N_GV_Top, prop_top, sep="__"), 
         diff_Bot=paste(N_GV_Bottom, prop_bottom, sep="__")) %>% 
  select(treatment, populationID, diff_Top, diff_Bot) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  ggplot(aes(group, Porportion, fill=group))+
  geom_boxplot()+
  ylim(value=c(0,0.35)) +
  scale_fill_manual(values=c("red3", "cornflowerblue"))+
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"))#+
  ggsave("plots/27_dominantFraction_boxplot_pink.png", width = 3, height = 4, dpi = 400)


data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_pool_top>(-0.15) & diff_pool_top<0.15)),
  #           N_GV_Pool=sum(diff_pool_top>0.15), 
  #           N_GV_Top=sum(diff_pool_top<(-0.15))*(-1),
  #           prop_Pool=N_GV_Pool/N_GV, 
  #           prop_Top=((N_GV_Top/N_GV))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_top>(-0.20) & diff_pool_top<0.20)),
            N_GV_Pool=sum(diff_pool_top>0.20), 
            N_GV_Top=sum(diff_pool_top<(-0.20))*(-1),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Top=((N_GV_Top/N_GV))) %>% 
  mutate(value_Pool=paste(N_GV_Pool, prop_Pool, sep="__"), 
         value_Top=paste(N_GV_Top, prop_Top, sep="__")) %>% 
  select(treatment, populationID, value_Pool, value_Top) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  # filter(group!="N_GV") %>% 
  # filter(group!="N_GV_NC") %>% 
  # filter(group %in% c("prop_top", "prop_bottom")) %>%
  # filter(group %in% c("N_GV_Pool", "N_GV_Top")) %>%
  ggplot(aes(populationID, Porportion, fill=group))+
  # geom_bar(stat="identity",position=position_dodge())+
  geom_bar(stat="identity")+
  geom_text(aes(label = abs(N_GV), 
                vjust = ifelse(Porportion >= 0, (-0.5), 1))) +
  scale_fill_manual(values=c("orange", "cornflowerblue"))+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  xlab("Population") +
  # ggplot(aes(group, NGenVar, fill=group))+
  # geom_boxplot()+
  # geom_jitter(alpha=0.9) +
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=90,size=12,colour="black",hjust= 1,vjust=0.5),
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/27_dominantFraction_pool_top_pink20.png", width = 7, height = 5, dpi = 400)


data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  # summarise(N_GV=n(), 
  #           N_GV_NC=sum((diff_pool_bottom>(-0.15) & diff_pool_bottom<0.15)),
  #           N_GV_Pool=sum(diff_pool_bottom>0.15), 
  #           N_GV_Bottom=sum(diff_pool_bottom<(-0.15))*(-1),
  #           prop_Pool=N_GV_Pool/N_GV, 
  #           prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_bottom>(-0.20) & diff_pool_bottom<0.20)),
            N_GV_Pool=sum(diff_pool_bottom>0.20), 
            N_GV_Bottom=sum(diff_pool_bottom<(-0.20))*(-1),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(value_Pool=paste(N_GV_Pool, prop_Pool, sep="__"), 
         value_Bottom=paste(N_GV_Bottom, prop_Bottom, sep="__")) %>% 
  select(treatment, populationID, value_Pool, value_Bottom) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  # filter(group!="N_GV") %>% 
  # filter(group!="N_GV_NC") %>% 
  # filter(group %in% c("prop_top", "prop_bottom")) %>%
  # filter(group %in% c("N_GV_Pool", "N_GV_Bottom")) %>%
  ggplot(aes(populationID, Porportion, fill=group))+
  # geom_bar(stat="identity",position=position_dodge())+
  geom_bar(stat="identity")+
  geom_text(aes(label = abs(N_GV), 
                vjust = ifelse(Porportion >= 0, (-0.5), 1))) +
  scale_fill_manual(values=c("red3", "orange"))+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  xlab("Population") +
  # ggplot(aes(group, NGenVar, fill=group))+
  # geom_boxplot()+
  # geom_jitter(alpha=0.9) +
  facet_grid(. ~ treatment, space="free", scale="free")+
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=90,size=12,colour="black",hjust= 1,vjust=0.5),
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/27_dominantFraction_pool_bottom_pink20.png", width = 7, height = 5, dpi = 400)



data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_top>(-0.20) & diff_pool_top<0.20)),
            N_GV_Pool=sum(diff_pool_top>0.20), 
            N_GV_Top=sum(diff_pool_top<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Top=((N_GV_Top/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Top=paste(N_GV_Top, prop_Top, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Top) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  ggplot(aes(group, Porportion, fill=group))+
  geom_boxplot()+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  scale_fill_manual(values=c("orange", "cornflowerblue"))+
  facet_grid(. ~ treatment, space="free", scale="free")+
  ylim(value=c(0,0.35)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave("plots/27_dominantFraction_boxplot_pool_top_pink.png", width = 3, height = 4, dpi = 400)


data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_bottom>(-0.20) & diff_pool_bottom<0.20)),
            N_GV_Pool=sum(diff_pool_bottom>0.20), 
            N_GV_Bottom=sum(diff_pool_bottom<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Bot=paste(N_GV_Bottom, prop_Bottom, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Bot) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  mutate(group=factor(group, levels=c("diff_Pool", "diff_Bot"))) %>% 
  ggplot(aes(group, Porportion, fill=group))+
  geom_boxplot()+
  # scale_fill_manual(values=c("red3", "cornflowerblue"))+
  scale_fill_manual(values=c("orange", "red3"))+
  facet_grid(. ~ treatment, space="free", scale="free")+
  ylim(value=c(0,0.35)) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("plots/27_dominantFraction_boxplot_pool_bottom_pink.png", width = 3, height = 4, dpi = 400)


# Statistical test:

# pink - top - bottom comparison
table_dominanFraction_top_bottom_pink<-data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>%
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_top_bottom>(-0.15) & diff_top_bottom<0.15)),
            N_GV_Top=sum(diff_top_bottom>0.15), 
            N_GV_Bottom=sum(diff_top_bottom<(-0.15)),
            prop_top=N_GV_Top/N_GV, 
            prop_bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Top=paste(N_GV_Top, prop_top, sep="__"), 
         diff_Bot=paste(N_GV_Bottom, prop_bottom, sep="__")) %>% 
  select(treatment, populationID, diff_Top, diff_Bot) %>% 
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  select(-populationID) %>%
  data.frame()

mod <- glm(Porportion~treatment/group,family=quasibinomial,table_dominanFraction_top_bottom_pink)
summary(mod)

# pink - Pool - top bottom comparison
table_dominanFraction_pool_top_pink<-data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_top>(-0.20) & diff_pool_top<0.20)),
            N_GV_Pool=sum(diff_pool_top>0.20), 
            N_GV_Top=sum(diff_pool_top<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Top=((N_GV_Top/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Top=paste(N_GV_Top, prop_Top, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Top) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  data.frame()

mod <- glm(Porportion~treatment/group,family=quasibinomial,table_dominanFraction_pool_top_pink)
summary(mod) 


# pink - Pool - Bottom bottom comparison
table_dominanFraction_pool_bottom_pink<-data.frame(treatment, populationID, Pool, Top, Bottom) %>% 
  filter(Pool!=0) %>%
  mutate(diff_pool_top=Pool-Top, 
         diff_pool_bottom=Pool-Bottom, 
         diff_top_bottom=Top-Bottom) %>%
  group_by(treatment, populationID) %>% 
  summarise(N_GV=n(), 
            N_GV_NC=sum((diff_pool_bottom>(-0.20) & diff_pool_bottom<0.20)),
            N_GV_Pool=sum(diff_pool_bottom>0.20), 
            N_GV_Bottom=sum(diff_pool_bottom<(-0.20)),
            prop_Pool=N_GV_Pool/N_GV, 
            prop_Bottom=((N_GV_Bottom/N_GV))) %>% 
  mutate(diff_Pool=paste(N_GV_Pool, prop_Pool, sep="__"),
         diff_Bot=paste(N_GV_Bottom, prop_Bottom, sep="__")) %>%
  select(treatment, populationID, diff_Pool, diff_Bot) %>%
  gather("group", "NGenVar", 3:4) %>% 
  separate(NGenVar, c("N_GV", "Porportion"), sep="__") %>% 
  mutate(Porportion=as.numeric(Porportion), 
         N_GV=as.numeric(N_GV)) %>% 
  mutate(group=factor(group, levels=c("diff_Pool", "diff_Bot"))) %>% 
  data.frame()

mod <- glm(Porportion~treatment/group,family=quasibinomial,table_dominanFraction_pool_bottom_pink)
summary(mod) 





#####
# Annotate GVs:
# vcf file to identify genomic effect: 

# red GVs
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))



table_vcf_allGV_snps_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(type_GV=="snps") %>% mutate(CHROM=Chrom, POS=Position, ID=".", 
                                     REF=Ref, ALT=VarAllele,  
                                     QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()


table_vcf_allGV_ins_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(type_GV=="indels" & grepl("\\+",VarAllele)) %>% 
  mutate(new_GV_code=paste0(Ref,gsub("\\+", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()

table_vcf_allGV_del_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(type_GV=="indels" & grepl("\\-",VarAllele)) %>% 
  mutate(new_GV_code=Ref, Ref=paste0(Ref,gsub("\\-", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()


table_vcf_allGV_formatted_red<-rbind(table_vcf_allGV_snps_red, 
                                     table_vcf_allGV_ins_red, 
                                     table_vcf_allGV_del_red) %>% 
  unique()

table_vcf_allGV_formatted_red[order(table_vcf_allGV_formatted_red$CHROM, 
                                    table_vcf_allGV_formatted_red$POS),] %>% 
  select(-GV_code) %>% 
  unique() %>% 
  write.table("./snpeff/table_vcf_allGV_formatted_red.vcf", quote = F, sep = "\t", row.names = FALSE)

##
# I did additional annotation file including only variants in parental phase III (time = 0) and at the end of the experiment (time = 53)

# red GVs
table_vcf_allGV_snps_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(time==0) %>% 
  filter(type_GV=="snps") %>% mutate(CHROM=Chrom, POS=Position, ID=".", 
                                     REF=Ref, ALT=VarAllele,  
                                     QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_ins_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(time==0) %>% 
  filter(type_GV=="indels" & grepl("\\+",VarAllele)) %>% 
  mutate(new_GV_code=paste0(Ref,gsub("\\+", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_del_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(time==0) %>% 
  filter(type_GV=="indels" & grepl("\\-",VarAllele)) %>% 
  mutate(new_GV_code=Ref, Ref=paste0(Ref,gsub("\\-", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_time0_formatted_red<-rbind(table_vcf_allGV_snps_red, 
                                           table_vcf_allGV_ins_red, 
                                           table_vcf_allGV_del_red) %>% 
  unique()
table_vcf_allGV_time0_formatted_red[order(table_vcf_allGV_time0_formatted_red$CHROM, 
                                          table_vcf_allGV_time0_formatted_red$POS),] %>% 
  select(-GV_code) %>% 
  unique() %>% 
  write.table("./snpeff/table_vcf_allGV_time0_formatted_red.vcf", quote = F, sep = "\t", row.names = FALSE)

# red GVs
table_vcf_allGV_snps_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(time==53) %>% 
  filter(time_groupGV!="phaseIII" | VarFreq>30) %>% 
  filter(type_GV=="snps") %>% mutate(CHROM=Chrom, POS=Position, ID=".", 
                                     REF=Ref, ALT=VarAllele,  
                                     QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_ins_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(time==53) %>% 
  filter(time_groupGV!="phaseIII" | VarFreq>30) %>% 
  filter(type_GV=="indels" & grepl("\\+",VarAllele)) %>% 
  mutate(new_GV_code=paste0(Ref,gsub("\\+", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_del_red<-phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3R") %>% 
  filter(time==53) %>% 
  filter(time_groupGV!="phaseIII" | VarFreq>30) %>% 
  filter(type_GV=="indels" & grepl("\\-",VarAllele)) %>% 
  mutate(new_GV_code=Ref, Ref=paste0(Ref,gsub("\\-", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_time53_formatted_red<-rbind(table_vcf_allGV_snps_red, 
                                            table_vcf_allGV_ins_red, 
                                            table_vcf_allGV_del_red) %>% 
  unique()
table_vcf_allGV_time53_formatted_red[order(table_vcf_allGV_time53_formatted_red$CHROM, 
                                           table_vcf_allGV_time53_formatted_red$POS),] %>% 
  select(-GV_code) %>% 
  unique() %>% 
  write.table("./snpeff/table_vcf_allGV_time53_formatted_red.vcf", quote = F, sep = "\t", row.names = FALSE)

# pink GVs
phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))



table_vcf_allGV_snps_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(type_GV=="snps") %>% mutate(CHROM=Chrom, POS=Position, ID=".", 
                                     REF=Ref, ALT=VarAllele,  
                                     QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()


table_vcf_allGV_ins_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(type_GV=="indels" & grepl("\\+",VarAllele)) %>% 
  mutate(new_GV_code=paste0(Ref,gsub("\\+", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()

table_vcf_allGV_del_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(type_GV=="indels" & grepl("\\-",VarAllele)) %>% 
  mutate(new_GV_code=Ref, Ref=paste0(Ref,gsub("\\-", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()


table_vcf_allGV_formatted_pink<-rbind(table_vcf_allGV_snps_pink, 
                                      table_vcf_allGV_ins_pink, 
                                      table_vcf_allGV_del_pink) %>% 
  unique()

table_vcf_allGV_formatted_pink[order(table_vcf_allGV_formatted_pink$CHROM, 
                                     table_vcf_allGV_formatted_pink$POS),] %>% 
  select(-GV_code) %>% 
  unique() %>% 
  write.table("./snpeff/table_vcf_allGV_formatted_pink.vcf", quote = F, sep = "\t", row.names = FALSE)

##
# I did additional annotation file including only variants in parental phase III (time = 0) and at the end of the experiment (time = 53)

# pink GVs
table_vcf_allGV_snps_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(time==0) %>% 
  filter(type_GV=="snps") %>% mutate(CHROM=Chrom, POS=Position, ID=".", 
                                     REF=Ref, ALT=VarAllele,  
                                     QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_ins_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(time==0) %>% 
  filter(type_GV=="indels" & grepl("\\+",VarAllele)) %>% 
  mutate(new_GV_code=paste0(Ref,gsub("\\+", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_del_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(time==0) %>% 
  filter(type_GV=="indels" & grepl("\\-",VarAllele)) %>% 
  mutate(new_GV_code=Ref, Ref=paste0(Ref,gsub("\\-", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_time0_formatted_pink<-rbind(table_vcf_allGV_snps_pink, 
                                            table_vcf_allGV_ins_pink, 
                                            table_vcf_allGV_del_pink) %>% 
  unique()
table_vcf_allGV_time0_formatted_pink[order(table_vcf_allGV_time0_formatted_pink$CHROM, 
                                           table_vcf_allGV_time0_formatted_pink$POS),] %>% 
  select(-GV_code) %>% 
  unique() %>% 
  write.table("./snpeff/table_vcf_allGV_time0_formatted_pink.vcf", quote = F, sep = "\t", row.names = FALSE)

# pink GVs
table_vcf_allGV_snps_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(time==53) %>% 
  filter(time_groupGV!="phaseIII" | VarFreq>30) %>% 
  filter(type_GV=="snps") %>% mutate(CHROM=Chrom, POS=Position, ID=".", 
                                     REF=Ref, ALT=VarAllele,  
                                     QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_ins_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(time==53) %>% 
  filter(time_groupGV!="phaseIII" | VarFreq>30) %>% 
  filter(type_GV=="indels" & grepl("\\+",VarAllele)) %>% 
  mutate(new_GV_code=paste0(Ref,gsub("\\+", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_del_pink<-phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) %>% 
  filter(parental=="P3P") %>% 
  filter(time==53) %>% 
  filter(time_groupGV!="phaseIII" | VarFreq>30) %>% 
  filter(type_GV=="indels" & grepl("\\-",VarAllele)) %>% 
  mutate(new_GV_code=Ref, Ref=paste0(Ref,gsub("\\-", "",VarAllele))) %>% 
  mutate(CHROM=Chrom, POS=Position, ID=".", 
         REF=Ref, ALT=new_GV_code, 
         QUAL=".", FILTER="PASS") %>% 
  select(GV_code, CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>% 
  unique()
table_vcf_allGV_time53_formatted_pink<-rbind(table_vcf_allGV_snps_pink, 
                                             table_vcf_allGV_ins_pink, 
                                             table_vcf_allGV_del_pink) %>% 
  unique()
table_vcf_allGV_time53_formatted_pink[order(table_vcf_allGV_time53_formatted_pink$CHROM, 
                                            table_vcf_allGV_time53_formatted_pink$POS),] %>% 
  select(-GV_code) %>% 
  unique() %>% 
  write.table("./snpeff/table_vcf_allGV_time53_formatted_pink.vcf", quote = F, sep = "\t", row.names = FALSE)

# vcf file were used in batch to run snpEff and produce eddited tables:

ann_all_GV_red<-read.table("./snpeff/ann_final_table_vcf_allGV_formatted_all_red_ed.txt", T) %>% mutate(GV_code=paste0(CHROM, "_", POS, "_", ALT)) %>% select(GV_code, Effect, Level, Gene)
ann_all_GV_pink<-read.table("./snpeff/ann_final_table_vcf_allGV_formatted_all_pink_ed.txt", T) %>% mutate(GV_code=paste0(CHROM, "_", POS, "_", ALT)) %>% select(GV_code, Effect, Level, Gene)


## 22_5 plots:
# distribution of effect size

# red

ann_all_GV_red<-read.table("./snpeff/ann_final_table_vcf_allGV_formatted_all_red_ed.txt", T) 

str(ann_all_GV_red)
ann_all_GV_snps_red <- ann_all_GV_red %>% 
  mutate(ref_char=nchar(as.vector(REF)), 
         alt_char=nchar(as.vector(ALT))) %>% 
  filter(ref_char==1, alt_char==1) %>% 
  mutate(GV_code=paste0(CHROM, "_", POS, "_", ALT)) %>% 
  select(GV_code, Effect, Level, Gene)

ann_all_GV_ins_red <- ann_all_GV_red %>% 
  mutate(ref_char=nchar(as.vector(REF)), 
         alt_char=nchar(as.vector(ALT))) %>% 
  filter(ref_char==1, alt_char>1) %>% 
  mutate(GV_code=paste0(CHROM, "_", POS, "_+",substring(ALT, ref_char+1))) %>% 
  select(GV_code, Effect, Level, Gene) 

ann_all_GV_del_red <- 
  ann_all_GV_red %>% 
  mutate(ref_char=nchar(as.vector(REF)), 
         alt_char=nchar(as.vector(ALT))) %>% 
  filter(ref_char>1, alt_char,ref_char) %>% 
  mutate(GV_code=paste0(CHROM, "_", POS, "_-",substring(REF, alt_char+1))) %>% 
  select(GV_code, Effect, Level, Gene) 

ann_all_GV_ed_red<-rbind(ann_all_GV_snps_red, ann_all_GV_ins_red, ann_all_GV_del_red)

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) %>% 
  merge(ann_all_GV_ed_red, by="GV_code", all.x=TRUE)

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_red <- phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) 

# write.table(phaseIII_table_total_red, "table_variants_red_18062020.txt", sep="\t", row.names = F)


phaseIII_table_total_red %>% 
  filter(time==0) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  ggplot(aes(VarFreq_max, fill=Level)) + 
  geom_histogram(position="dodge") +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), #element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  facet_grid(Effect ~ .) #+  
  ggsave(paste0("plots/22_5_effect_parental_phaseIII_red.png"), width = 6, height = 9, dpi = 400)


phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  ggplot(aes(VarFreq_max, fill=Level)) + 
  geom_histogram(position="dodge") +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), #element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=10)) +
  facet_grid(Effect ~ ., scale="free") #+  
  ggsave(paste0("plots/22_5_effect_phaseIII_red.png"), width = 10, height = 15, dpi = 400)


phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  # filter(time_groupGV!="phaseIII") %>%
  ungroup() %>%
  mutate(VarFreq_max=VarFreq_max/100) %>%
  mutate(Level=factor(Level, levels=c("LOW", "MODERATE", "MODIFIER", "HIGH"))) %>%
  ggplot(aes(Level, VarFreq_max, fill=Level)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.3) +
  scale_fill_manual(values=rev(heat.colors(5))) +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=10)) #+
  # ggsave(paste0("plots/22_5_effect_phaseIII_boxplot_withoutNewPhaseIII_red.png"), width = 5, height = 6, dpi = 400)
  ggsave(paste0("plots/22_5_effect_phaseIII_boxplot_red.png"), width = 5, height = 6, dpi = 400)

# test: 
test_table_red<-phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  # filter(time_groupGV!="phaseIII") %>%
  ungroup() %>%
  mutate(VarFreq_max=VarFreq_max/100) %>%
  mutate(Level=factor(Level, levels=c("LOW", "MODERATE", "MODIFIER", "HIGH"))) %>%
  select(Level, VarFreq_max) %>%
  data.frame()

head(test_table_red)
mod <- glm(VarFreq_max~Level,family=binomial(link=logit),test_table_red)
summary(mod)
mod <- glm(VarFreq_max~Level,family=quasibinomial(link=logit),test_table_red)
summary(mod)


phaseIII_table_total_red %>% 
  filter(time %in% c(0,53)) %>% 
  #mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time,pop_treatment,treatment, VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time,pop_treatment,treatment, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  select(Chrom,Position,GV_code,time_groupGV, time,pop_treatment,treatment,VarFreq, Effect, Level, Gene) %>% 
  # spread(time_pop_ID, VarFreq, fill=0) %>% 
  # gather("time_pop_ID", "VarFreq", 8:73) %>% 
  #tidyr::separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% head()
  filter(time_groupGV!="phaseIII") %>% 
  group_by(GV_code, Effect, Level, Gene, treatment) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  mutate(VarFreq_max=VarFreq_max/100) %>%
  ungroup() %>% 
  mutate(Level=factor(Level, levels=c("LOW", "MODERATE", "MODIFIER", "HIGH"))) %>%
  mutate(treatment=factor(treatment, levels=c("P3R", "Allo_T", "Allo_B", 
                                              "Para_T", "Para_B", 
                                              "LM", "SYM"))) %>%
  ggplot(aes(Level, VarFreq_max, fill=Level)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.3) +
  scale_fill_manual(values=rev(heat.colors(5))) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  facet_grid(. ~ treatment, scale="free") +
  # ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_red.png"), width = 11, height = 5, dpi = 400)
  # ggsave(paste0("22_5_effect_phaseIII_bytreatment_boxplot_withoutPhaseIII_red.png"), width = 11, height = 5, dpi = 400)

phaseIII_table_total_red %>% 
  filter(time %in% c(0,53)) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  group_by(GV_code, Effect, Level, Gene, treatment) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  ggplot(aes(Level, VarFreq_max, fill=Level)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.3) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  facet_grid(. ~ treatment, scale="free") #+
  # ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_withoutPhaseIII_red.png"), width = 11, height = 5, dpi = 400)
  ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_red.png"), width = 11, height = 5, dpi = 400)

#pink

ann_all_GV_pink<-read.table("./snpeff/ann_final_table_vcf_allGV_formatted_all_pink_ed.txt", T) 

ann_all_GV_snps_pink <- ann_all_GV_pink %>% 
  mutate(ref_char=nchar(as.vector(REF)), 
         alt_char=nchar(as.vector(ALT))) %>% 
  filter(ref_char==1, alt_char==1) %>% 
  mutate(GV_code=paste0(CHROM, "_", POS, "_", ALT)) %>% 
  select(GV_code, Effect, Level, Gene)

ann_all_GV_ins_pink <- ann_all_GV_pink %>% 
  mutate(ref_char=nchar(as.vector(REF)), 
         alt_char=nchar(as.vector(ALT))) %>% 
  filter(ref_char==1, alt_char>1) %>% 
  mutate(GV_code=paste0(CHROM, "_", POS, "_+",substring(ALT, ref_char+1))) %>% 
  select(GV_code, Effect, Level, Gene) 

ann_all_GV_del_pink <- 
  ann_all_GV_pink %>% 
  mutate(ref_char=nchar(as.vector(REF)), 
         alt_char=nchar(as.vector(ALT))) %>% 
  filter(ref_char>1, alt_char,ref_char) %>% 
  mutate(GV_code=paste0(CHROM, "_", POS, "_-",substring(REF, alt_char+1))) %>% 
  select(GV_code, Effect, Level, Gene) 

ann_all_GV_ed_pink<-rbind(ann_all_GV_snps_pink, ann_all_GV_ins_pink, ann_all_GV_del_pink)

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) %>% 
  merge(ann_all_GV_ed_pink, by="GV_code", all.x=TRUE)

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_pink <- phaseIII_table_total_pink %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) 

# write.table(phaseIII_table_total_red, "table_variants_red_18062020.txt", sep="\t", 
#             row.names = F, quote = F)
# write.table(phaseIII_table_total_pink, "table_variants_pink_18062020.txt", sep="\t", 
#             row.names = F, quote = F)

phaseIII_table_total_pink %>% 
  filter(time==0) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  ggplot(aes(VarFreq_max, fill=Level)) + 
  geom_histogram(position="dodge") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), #element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=10)) +
  facet_grid(Effect ~ .) #+
  # ggsave(paste0("plots/22_5_effect_parental_phaseIII_pink.png"), width = 7, height = 9, dpi = 400)


phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  ggplot(aes(VarFreq_max, fill=Level)) + 
  geom_histogram(position="dodge") +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), #element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=10)) +
  facet_grid(Effect ~ ., scale="free") #+  
  # ggsave(paste0("plots/22_5_effect_phaseIII_pink.png"), width = 10, height = 15, dpi = 400)


phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  filter(!(is.na(Level))) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  # filter(time_groupGV!="phaseIII") %>%
  ungroup() %>%
  mutate(VarFreq_max=VarFreq_max/100) %>%
  mutate(Level=factor(Level, levels=c("LOW", "MODERATE", "MODIFIER", "HIGH"))) %>%
  ggplot(aes(Level, VarFreq_max, fill=Level)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.3) +
  scale_fill_manual(values=rev(heat.colors(5))) +
  #geom_histogram(position="dodge")
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=10)) #+
  # ggsave(paste0("plots/22_5_effect_phaseIII_boxplot_withoutNewPhaseIII_pink.png"), width = 5, height = 6, dpi = 400)
  # ggsave(paste0("plots/22_5_effect_phaseIII_boxplot_pink.png"), width = 5, height = 6, dpi = 400)


# test: 
test_table_pink<-phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  filter(!(is.na(Level))) %>% 
  group_by(GV_code, Effect, Level, Gene, time_groupGV) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  filter(time_groupGV!="phaseIII") %>%
  ungroup() %>%
  mutate(VarFreq_max=VarFreq_max/100) %>%
  mutate(Level=factor(Level, levels=c("LOW", "MODERATE", "MODIFIER", "HIGH"))) %>%
  select(Level, VarFreq_max) %>%
  data.frame()

head(test_table_pink)
mod <- glm(VarFreq_max~Level,family=binomial(link=logit),test_table_pink)
summary(mod)
mod <- glm(VarFreq_max~Level,family=quasibinomial(link=logit),test_table_pink)
summary(mod)


phaseIII_table_total_pink %>% 
  filter(time %in% c(0,53)) %>% 
  filter(!(is.na(Level))) %>% 
  #mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time,pop_treatment,treatment, VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time,pop_treatment,treatment, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom,Position,GV_code,time_groupGV, time,pop_treatment,treatment,VarFreq, Effect, Level, Gene) %>% 
  filter(time_groupGV!="phaseIII") %>%
  group_by(GV_code, Effect, Level, Gene, treatment) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  mutate(VarFreq_max=VarFreq_max/100) %>%
  ungroup() %>% 
  mutate(Level=factor(Level, levels=c("LOW", "MODERATE", "MODIFIER", "HIGH"))) %>%
  mutate(treatment=factor(treatment, levels=c("P3P", "Allo_T", "Allo_B", 
                                              "Para_T", "Para_B", 
                                              "LM", "SYM"))) %>%
  ggplot(aes(Level, VarFreq_max, fill=Level)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.3) +
  scale_fill_manual(values=rev(heat.colors(5))) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  facet_grid(. ~ treatment, scale="free") +
  # ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_pink.png"), width = 11, height = 5, dpi = 400)
  ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_withoutPhaseIII_pink.png"), width = 11, height = 5, dpi = 400)
  


phaseIII_table_total_pink %>% 
  filter(time %in% c(0,53)) %>% 
  filter(!(is.na(Level))) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Effect, Level, Gene) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 8:73) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  # filter(time_groupGV!="phaseIII") %>% 
  group_by(GV_code, Effect, Level, Gene, treatment) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  #group_by(Effect, time_groupGV) %>% 
  #dplyr::summarize(N_GVs=n(), VarFreq=mean(VarFreq)) %>% 
  ggplot(aes(VarFreq_max, fill=Level)) + 
  geom_vline(xintercept = 0) +
  geom_histogram(position="dodge") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), #element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  facet_grid(Effect ~ treatment, scale="free") #+  
  # ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_pink.png"), width = 11, height = 10, dpi = 400)


phaseIII_table_total_pink %>% 
  filter(time %in% c(0,53)) %>% 
  filter(!(is.na(Level))) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  group_by(GV_code, Effect, Level, Gene, treatment) %>% 
  dplyr::summarize(N_population=n(), VarFreq_max=max(VarFreq)) %>% 
  ggplot(aes(Level, VarFreq_max, fill=Level)) + 
  geom_boxplot() +
  geom_point(position="jitter", alpha=0.3) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  facet_grid(. ~ treatment, scale="free") +
  # ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_withoutPhaseIII_pink.png"), width = 11, height = 5, dpi = 400)
  ggsave(paste0("plots/22_5_effect_phaseIII_bytreatment_boxplot_pink.png"), width = 11, height = 5, dpi = 400)

#
##### 
# GVs per gene:

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) %>% 
  merge(ann_all_GV_ed_red, by="GV_code", all.x=TRUE)

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) %>% 
  merge(ann_all_GV_ed_pink, by="GV_code", all.x=TRUE)


phaseIII_table_total_both <- rbind(phaseIII_table_total_red, phaseIII_table_total_pink) 

phaseIII_table_total_both$time_groupGV<-"phaseIII"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
# using high fq variants in parentals phase I
# phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_parental_GV_red]<-"polymorphic_parental"
# phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_parental_GV_pink]<-"polymorphic_parental"
# phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% fix_parental_GV_red]<-"fix_parental"
# phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% fix_parental_GV_pink]<-"fix_parental"
# using min 20% variants parentals phase I
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_both$time_groupGV[phaseIII_table_total_both$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"

phaseIII_table_total_both$time_groupGV <- factor(phaseIII_table_total_both$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))


phaseIII_table_total_both %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polimorphic_parental"))) %>% 
  filter(time==53) %>% 
  unique() %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, treatment, pop_treatment, 
           Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, Effect, Level, Gene) %>% 
  dplyr::summarize(N_populations=n(), VarFreq_max=max(VarFreq)) %>% 
  filter(!(Effect %in% c("intergenic_region", "3_prime_UTR_variant", "5_prime_UTR_variant"))) %>% 
  group_by(Gene) %>% 
  filter(!(is.na(Gene))) %>% 
  #filter(Gene!="SPBPJ4664.02") %>% 
  #filter(Gene=="mal1") -> mal1
  dplyr::summarize(N_GVs=n()) %>% 
  # filter(N_GVs>10)
  # filter(N_GVs>1) %>% 
  # arrange(N_GVs) %>%
  # write.table(file = "text_files/table_Num_GVs_perGene_both.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  filter(N_GVs<50) %>%
  ggplot(aes(N_GVs)) +
  geom_histogram(fill="black", binwidth = 1) + 
  # xlim(c(0,20)) +
  theme_classic() + 
  theme(panel.border = element_blank(),  
        panel.grid.major = element_blank(), #element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
ggsave(paste0("plots/22_6_distributionGVs_perGene_phaseIII_both.png"), 
       width = 5, height = 4, dpi = 400)



phaseIII_table_total_both %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polimorphic_parental"))) %>% 
  filter(time==53) %>% 
  unique() %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, treatment, pop_treatment, 
           Effect, Level, Gene) %>% 
  #filter(GV_code=="II_4376054_G") %>% data.frame()
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, Effect, Level, Gene) %>% 
  dplyr::summarize(N_populations=n(), VarFreq_max=max(VarFreq)) %>% 
  # filter(Gene=="SPBPJ4664.02") %>%
  # filter(Gene=="mal1") %>%
  #filter(Gene=="pfl4") 
  #filter(Gene=="map4") 
  filter(Gene=="pfl2") %>% 
  #filter(Gene=="SPCC1235.01") 
  #filter(Gene=="ubr11") 
  #filter(Gene=="ght8") 
  #filter(Gene=="pvg3") #*
  #filter(Gene=="pvg2") #*
  #filter(Gene=="pvg1") #*
#filter(Gene=="pvg5") #*
#filter(Gene=="efc25")
#filter(Gene=="msa1")
#filter(Gene=="meu7") 
#filter(Gene=="agn1") 
#filter(Gene=="SPAPB2C8.01")
#filter(Gene=="ubi4")
#filter(Gene=="sck1") 
# filter(Gene=="wtf21")
ggplot(aes(N_populations, fill=time_groupGV)) +
  geom_histogram(position="dodge") + 
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  ggsave(paste0("plots/22_6_distributionGVs_pfl2_phaseIII_both.png"), width = 8, height = 4, dpi = 400)

phaseIII_table_total_both %>% 
  filter(!(time_groupGV %in% c("fix_parental", "polimorphic_parental"))) %>% 
  filter(time==53) %>% 
  unique() %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, treatment, pop_treatment, 
           Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, Effect, Level, Gene) %>% 
  dplyr::summarize(N_populations=n(), VarFreq_max=max(VarFreq)) %>% 
  filter(!(Effect %in% c("intergenic_region", "3_prime_UTR_variant", "5_prime_UTR_variant"))) %>% 
  group_by(Gene) %>% 
  filter(!(is.na(Gene))) %>% 
  #filter(Gene!="SPBPJ4664.02") %>% 
  #filter(Gene=="mal1") -> mal1
  dplyr::summarize(N_GVs=n()) %>% 
  head()








#####
# Logistic Regression with binomial/quasibinomial model - Allopatric treatment

# red populations:
simple_mean_coverage_table<-sample_information %>% 
  filter(!(is.na(time))) %>% 
  mutate(time_pop_ID=paste0(time,"__",population, "_", treatment,"__",treatment)) %>% 
  select(time_pop_ID, total_mean_coverage) 

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) %>% 
  merge(ann_all_GV_ed_red, by="GV_code", all.x=TRUE)

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))
phaseIII_table_total_red <- phaseIII_table_total_red %>% 
  filter(!(time_groupGV %in% c("fix_parental"))) 
  # filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) 

phaseIII_table_total_red_test <- phaseIII_table_total_red %>%
  filter(time %in% c(53)) %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Effect!="intergenic_region") %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  mutate(VarFreq_ID=paste(Reads1+Reads2,VarFreq, sep="__")) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, 
         VarFreq_ID, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(VarFreq_ID=VarFreq_ID[which(Pvalue==min(Pvalue))]) %>% 
  spread(time_pop_ID, VarFreq_ID, fill="0__0") %>% 
  gather("time_pop_ID", "VarFreq_ID", 8:29) %>% 
  merge(simple_mean_coverage_table, all.x=T) %>% 
  mutate(VarFreq_ID_ed=ifelse(VarFreq_ID=="0__0", 
                              paste(total_mean_coverage, 0, sep="__"), 
                              VarFreq_ID)) %>% 
  select(-VarFreq_ID, -total_mean_coverage) 

Allo_phaseIII_table_total_red_test <- phaseIII_table_total_red_test %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  ungroup() %>% 
  select(GV_code, pop_ID, treatment, VarFreq_ID_ed) %>%
  spread(GV_code, VarFreq_ID_ed) %>% 
  ungroup() %>% 
  mutate(treatment=as.factor(treatment)) %>% 
  select(-pop_ID) 

GV_code<-c()
pvalue_binomial<-c()
pvalue_quasibinomial<-c()
Allo_B<-c()
Allo_T<-c()
for (i in seq(2,dim(Allo_phaseIII_table_total_red_test)[2])){
  print(i)
  GV_code[i-1]<-names(Allo_phaseIII_table_total_red_test[i])
  print(GV_code[i-1])
  #Allo_phaseIII_table_total_red_test
  variant_table<-Allo_phaseIII_table_total_red_test[,c(1, i)] %>% 
    separate(names(Allo_phaseIII_table_total_red_test)[i], c("cov", "fq"), sep="__") %>% 
    mutate(alt_counts=round(as.numeric(cov)*as.numeric(fq)/100),
           ref_counts=as.numeric(cov)-alt_counts) %>%
    select(treatment, alt_counts, ref_counts) 
  #print(head(variant_table))
  mod <- glm(cbind(alt_counts,ref_counts)~treatment,family=binomial(link=logit),variant_table)
  summury_mod<-summary(mod)
  Allo_B[i-1]<-as.vector(summury_mod$coefficients[1,1])
  Allo_T[i-1]<-as.vector(summury_mod$coefficients[2,1])
  pvalue_binomial[i-1]<-as.vector(summury_mod$coefficients[2,4])
  mod <- glm(cbind(alt_counts,ref_counts)~treatment,family=quasibinomial(link=logit),variant_table)
  summury_mod<-summary(mod)
  pvalue_quasibinomial[i-1]<-as.vector(summury_mod$coefficients[2,4])
}


glm_test_allo_red<-data.frame(GV_code, Allo_B, Allo_T, pvalue_binomial, cor_fdr_pval_binomial=p.adjust(pvalue_binomial, method = "fdr"), pvalue_quasibinomial, cor_fdr_pval_quasibinomial=p.adjust(pvalue_quasibinomial, method = "fdr"))

write.table(glm_test_allo_red, file = "text_files/table_glm_test_allo_red.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# pink populations:

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) %>% 
  merge(ann_all_GV_ed_pink, by="GV_code", all.x=TRUE)

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

phaseIII_table_total_pink <- phaseIII_table_total_pink %>%
  filter(!(time_groupGV %in% c("fix_parental")))
  # filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental")))

phaseIII_table_total_pink_test <- phaseIII_table_total_pink %>%
  filter(time %in% c(53)) %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Effect!="intergenic_region") %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  mutate(VarFreq_ID=paste(Reads1+Reads2,VarFreq, sep="__")) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, 
         VarFreq_ID, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(VarFreq_ID=VarFreq_ID[which(Pvalue==min(Pvalue))[1]]) %>% 
  spread(time_pop_ID, VarFreq_ID, fill="0__0") %>% 
  gather("time_pop_ID", "VarFreq_ID", 8:28) %>% 
  merge(simple_mean_coverage_table, all.x=T) %>% 
  mutate(VarFreq_ID_ed=ifelse(VarFreq_ID=="0__0", 
                              paste(total_mean_coverage, 0, sep="__"), 
                              VarFreq_ID)) %>% 
  select(-VarFreq_ID, -total_mean_coverage) 

Allo_phaseIII_table_total_pink_test <- phaseIII_table_total_pink_test %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  ungroup() %>% 
  select(GV_code, pop_ID, treatment, VarFreq_ID_ed) %>%
  spread(GV_code, VarFreq_ID_ed) %>% 
  ungroup() %>% 
  mutate(treatment=as.factor(treatment)) %>% 
  select(-pop_ID) 

GV_code<-c()
pvalue_binomial<-c()
pvalue_quasibinomial<-c()
Allo_B<-c()
Allo_T<-c()
for (i in seq(2,dim(Allo_phaseIII_table_total_pink_test)[2])){
  print(i)
  GV_code[i-1]<-names(Allo_phaseIII_table_total_pink_test[i])
  print(GV_code[i-1])
  #Allo_phaseIII_table_total_pink_test
  variant_table<-Allo_phaseIII_table_total_pink_test[,c(1, i)] %>% 
    separate(names(Allo_phaseIII_table_total_pink_test)[i], c("cov", "fq"), sep="__") %>% 
    mutate(alt_counts=round(as.numeric(cov)*as.numeric(fq)/100),
           ref_counts=as.numeric(cov)-alt_counts) %>%
    select(treatment, alt_counts, ref_counts) 
  #print(head(variant_table))
  mod <- glm(cbind(alt_counts,ref_counts)~treatment,family=binomial(link=logit),variant_table)
  summury_mod<-summary(mod)
  Allo_B[i-1]<-as.vector(summury_mod$coefficients[1,1])
  Allo_T[i-1]<-as.vector(summury_mod$coefficients[2,1])
  pvalue_binomial[i-1]<-as.vector(summury_mod$coefficients[2,4])
  mod <- glm(cbind(alt_counts,ref_counts)~treatment,family=quasibinomial(link=logit),variant_table)
  summury_mod<-summary(mod)
  pvalue_quasibinomial[i-1]<-as.vector(summury_mod$coefficients[2,4])
}


glm_test_allo_pink<-data.frame(GV_code, Allo_B, Allo_T, pvalue_binomial, cor_fdr_pval_binomial=p.adjust(pvalue_binomial, method = "fdr"), pvalue_quasibinomial, cor_fdr_pval_quasibinomial=p.adjust(pvalue_quasibinomial, method = "fdr"))

write.table(glm_test_allo_pink, file = "text_files/table_glm_test_allo_pink.txt", quote = FALSE, sep = "\t", row.names = FALSE)

top_bottom_diffGV_pink<-glm_test_allo_pink %>% 
  filter(cor_fdr_pval_binomial<0.05 | pvalue_quasibinomial<0.05) %>% 
  select(GV_code)

top_bottom_diffGV_red<-glm_test_allo_red %>% 
  filter(cor_fdr_pval_binomial<0.05 | pvalue_quasibinomial<0.05) %>% 
  select(GV_code)

glm_test_allo_pink %>% 
  # filter(pvalue_binomial<0.1) %>% 
  # filter(pvalue_quasibinomial<0.1) %>% 
  filter(pvalue_binomial<0.05) %>% 
  ggplot(aes(cor_fdr_pval_binomial)) +
  geom_histogram()
# head()
filter(pvalue_binomial<0.05) %>% 
  # filter(pvalue_quasibinomial<0.05) %>% 
  # filter(cor_fdr_pval_quasibinomial<0.05)
  # filter(pvalue_quasibinomial<0.05) %>% 
  select(GV_code)


dim(glm_test_allo_pink)

glm_test_allo_red %>% 
  #filter(cor_fdr_pval_quasibinomial<0.05)
  # filter(pvalue_binomial<0.05) 
  filter(pvalue_quasibinomial<0.05) %>% 
  select(GV_code)

glm_test_allo_red %>% 
  filter(pvalue_quasibinomial<0.05) 

glm_test_allo_pink %>% 
  filter(pvalue_quasibinomial<0.05) 


#####
# Changes in frequency - phase III
# Allopatric populations:
# Red
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) %>% 
  merge(ann_all_GV_ed_red, by="GV_code", all.x=TRUE)

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))
# phaseIII_table_total_red <- phaseIII_table_total_red %>% 
#  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) 

parental_phaseIII_table_total_red_top_bottom_diffGV <- phaseIII_table_total_red %>% 
  filter(GV_code %in% top_bottom_diffGV_red$GV_code) %>% 
  filter(time %in% c(0)) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq_parental=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  ungroup() %>% 
  select(GV_code, VarFreq_parental) 

Num_GV_perGene<-read.table("./text_files/table_Num_GVs_perGene_allgenes_both.txt", T) %>% mutate(N_GVs_gene=N_GVs) %>% select(-N_GVs)

# phaseIII_table_total_red_top_bottom_diffGV <- phaseIII_table_total_red %>% 
#   filter(GV_code %in% top_bottom_diffGV_red$GV_code)


phaseIII_table_total_red_top_bottom_diffGV <- 
  phaseIII_table_total_red %>% 
  filter(GV_code %in% top_bottom_diffGV_red$GV_code) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, 
         VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(-N_GV) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 8:72) %>% 
  merge(parental_phaseIII_table_total_red_top_bottom_diffGV, all=TRUE) %>% 
  merge(Num_GV_perGene, all.x = TRUE)

phaseIII_table_total_red_top_bottom_diffGV$VarFreq_parental[is.na(phaseIII_table_total_red_top_bottom_diffGV$VarFreq_parental)]<-0

phaseIII_table_total_private_GVs_allo_red_top_bottom_diffGV <-
  phaseIII_table_total_red_top_bottom_diffGV %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  # filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  group_by(GV_code) %>% 
  filter(VarFreq>30) %>% 
  dplyr::summarize(N_populations=n()) %>% 
  # ggplot(aes(N_populations)) +
  # geom_histogram() 
  # xlim(c(0,10))
  filter(N_populations>4) %>%
  ungroup() %>% 
  select(GV_code)
  

phaseIII_table_total_red_top_bottom_diffGV %>% 
  mutate(changeVarFreq=VarFreq-VarFreq_parental) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  filter(GV_code %in% phaseIII_table_total_private_GVs_allo_red_top_bottom_diffGV$GV_code) %>% 
  filter(GV_code %in% top_bottom_diffGV_red$GV_code) %>% 
  #filter(Chrom==chr) %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Position==3870142)
  #filter(Gene!="SPBPJ4664.02") %>%
  separate(GV_code, c("chr", "pos", "var"), "_") %>% 
  mutate(var=sub("(^....).*","\\1",var), GV_code=paste0(chr, ":", pos, "_", var)) %>% 
  mutate(GV_code=factor(GV_code, levels=c(unique(c(GV_code[sort.list(Position)]))))) %>% 
  ggplot(aes(GV_code, VarFreq/100, fill=treatment)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  #geom_point(alpha=0.3) +
  geom_point(aes(fill = treatment), size = 1.5, shape=21, alpha=0.7, position = position_jitterdodge()) +
  geom_point(aes(GV_code, VarFreq_parental/100), colour="orange", size=2, alpha=0.4) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  #scale_fill_brewer(palette="Set1") +
  #coord_flip() +
  ylab("Frequency Derived Allele") +
  xlab("Variant ID") +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="top", 
        legend.text = element_text(size=12, colour="black"),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title = element_text(size=13, colour="black"), 
        # strip.text.y = element_text(angle = 360, size=12), 
        strip.text = element_blank(), 
        strip.background = element_blank())  +
  ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_1_red.png"), width = 12, height = 4, dpi = 400)
  # ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_2_red.png"), width = 37, height = 7, dpi = 400)
  # ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_3_red.png"), width = 17, height = 7, dpi = 400)

# by gene:
  

phaseIII_table_total_bygene_red_top_bottom_diffGV<-
  phaseIII_table_total_red_top_bottom_diffGV %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  filter(time==53) %>% 
  filter(GV_code %in% phaseIII_table_total_private_GVs_allo_red_top_bottom_diffGV$GV_code) %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Effect!="intergenic_region") %>% 
  mutate(time_pop_ID=paste(time, pop_ID, treatment, sep="__")) %>% 
  ungroup() %>% 
  select(time_pop_ID, Chrom, Gene, VarFreq) %>% 
  group_by(time_pop_ID, Chrom, Gene) %>% 
  dplyr::summarize(VarFreq_max=max(VarFreq)) %>%
  merge((phaseIII_table_total_red_top_bottom_diffGV %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  filter(time==53) %>% 
  filter(GV_code %in% phaseIII_table_total_private_GVs_allo_red_top_bottom_diffGV$GV_code) %>% 
  group_by(Gene) %>% 
  dplyr::summarize(N_GVs=length(unique(Position)), Position_min=min(Position))), all.x = T) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") 

phaseIII_table_total_bygene_red_top_bottom_diffGV<-
  cbind(phaseIII_table_total_bygene_red_top_bottom_diffGV, 
        Gene_NumGV=factor(paste0(phaseIII_table_total_bygene_red_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$Position_min), 
                          levels=unique(paste0(phaseIII_table_total_bygene_red_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$Position_min)[order(phaseIII_table_total_bygene_red_top_bottom_diffGV$Position_min)])))

phaseIII_table_total_bygene_red_top_bottom_diffGV %>% 
  ggplot(aes(Gene_NumGV, VarFreq_max, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, shape=21, alpha=0.3, position = position_jitterdodge()) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
ggsave(paste0("plots/22_8_changeFq_Allo_phaseIII_byGene_red.png"), width = 15, height = 8, dpi = 400)



# pink
phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) %>% 
  merge(ann_all_GV_ed_pink, by="GV_code", all.x=TRUE)

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))
# phaseIII_table_total_pink <- phaseIII_table_total_pink %>% 
#  filter(!(time_groupGV %in% c("fix_parental", "polymorphic_parental"))) 

# phaseIII_table_total_pink_top_bottom_diffGV <- phaseIII_table_total_pink %>% 
#   filter(GV_code %in% top_bottom_diffGV_pink$GV_code)

parental_phaseIII_table_total_pink_top_bottom_diffGV <- phaseIII_table_total_pink %>% 
  filter(GV_code %in% top_bottom_diffGV_pink$GV_code) %>% 
  filter(time %in% c(0)) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV,time_pop_ID,VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq_parental=VarFreq[which(Pvalue==min(Pvalue))]) %>% 
  ungroup() %>% 
  select(GV_code, VarFreq_parental) 

Num_GV_perGene<-read.table("./text_files/table_Num_GVs_perGene_allgenes_both.txt", T) %>% mutate(N_GVs_gene=N_GVs) %>% select(-N_GVs)

phaseIII_table_total_pink_top_bottom_diffGV <-phaseIII_table_total_pink %>% 
  filter(GV_code %in% top_bottom_diffGV_pink$GV_code) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, 
         VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(-N_GV) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 8:72) %>% 
  merge(parental_phaseIII_table_total_pink_top_bottom_diffGV, all.x=TRUE) %>% 
  merge(Num_GV_perGene, all.x = TRUE)

phaseIII_table_total_pink_top_bottom_diffGV$VarFreq_parental[is.na(phaseIII_table_total_pink_top_bottom_diffGV$VarFreq_parental)]<-0

phaseIII_table_total_private_GVs_allo_pink_top_bottom_diffGV <-
  phaseIII_table_total_pink_top_bottom_diffGV %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  # filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  group_by(GV_code) %>% 
  filter(VarFreq>20) %>% 
  dplyr::summarize(N_populations=n()) %>% 
  # ggplot(aes(N_populations)) +
  # geom_histogram() 
  # xlim(c(0,10))
  filter(N_populations>4) %>% 
  ungroup() %>% 
  select(GV_code)


phaseIII_table_total_pink_top_bottom_diffGV %>% 
  mutate(changeVarFreq=VarFreq-VarFreq_parental) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  filter(GV_code %in% phaseIII_table_total_private_GVs_allo_pink_top_bottom_diffGV$GV_code) %>% 
  #filter(Chrom==chr) %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Position==3870142)
  #filter(Gene!="SPBPJ4664.02") %>%
  separate(GV_code, c("chr", "pos", "var"), "_") %>% 
  mutate(var=sub("(^....).*","\\1",var), GV_code=paste0(chr, ":", pos, "_", var)) %>% 
  mutate(GV_code=factor(GV_code, levels=c(unique(c(GV_code[sort.list(Position)]))))) %>% 
  ggplot(aes(GV_code, VarFreq/100, fill=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_point(alpha=0.3) +
  geom_point(aes(fill = treatment), size = 1, shape=21, alpha=0.5, position = position_jitterdodge()) +
  geom_point(aes(GV_code, VarFreq_parental/100), colour="orange", size=2, alpha=0.4) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  #scale_fill_brewer(palette="Set1") +
  #coord_flip() +
  ylab("Frequency Derived Allele") +
  xlab("Variant ID") +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="top", 
        legend.text = element_text(size=12, colour="black"),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title = element_text(size=13, colour="black"), 
        # strip.text.y = element_text(angle = 360, size=12), 
        strip.text = element_blank(), 
        strip.background = element_blank())  +
ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_1_pink.png"), width = 15, height = 5, dpi = 400)
# ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_2_pink.png"), width = 37, height = 7, dpi = 400)
# ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_3_pink.png"), width = 17, height = 7, dpi = 400)


phaseIII_table_total_red_top_bottom_diffGV %>% 
  mutate(changeVarFreq=VarFreq-VarFreq_parental) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  filter(GV_code %in% phaseIII_table_total_private_GVs_allo_red_top_bottom_diffGV$GV_code) %>% 
  filter(GV_code %in% top_bottom_diffGV_red$GV_code) %>% 
  #filter(Chrom==chr) %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Position==3870142)
  #filter(Gene!="SPBPJ4664.02") %>%
  separate(GV_code, c("chr", "pos", "var"), "_") %>% 
  mutate(var=sub("(^....).*","\\1",var), GV_code=paste0(chr, ":", pos, "_", var)) %>% 
  mutate(GV_code=factor(GV_code, levels=c(unique(c(GV_code[sort.list(Position)]))))) %>% 
  ggplot(aes(GV_code, VarFreq/100, fill=treatment)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  #geom_point(alpha=0.3) +
  geom_point(aes(fill = treatment), size = 1.5, shape=21, alpha=0.7, position = position_jitterdodge()) +
  geom_point(aes(GV_code, VarFreq_parental/100), colour="orange", size=2, alpha=0.4) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  #scale_fill_brewer(palette="Set1") +
  #coord_flip() +
  ylab("Frequency Derived Allele") +
  xlab("Variant ID") +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="top", 
        legend.text = element_text(size=12, colour="black"),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title = element_text(size=13, colour="black"), 
        # strip.text.y = element_text(angle = 360, size=12), 
        strip.text = element_blank(), 
        strip.background = element_blank())  +
  ggsave(paste0("plots/22_7_changeFq_Allo_phaseIII_1_red.png"), width = 12, height = 4, dpi = 400)

# by gene:

phaseIII_table_total_bygene_pink_top_bottom_diffGV<-
  phaseIII_table_total_pink_top_bottom_diffGV %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
  filter(time==53) %>% 
  filter(GV_code %in% phaseIII_table_total_private_GVs_allo_pink_top_bottom_diffGV$GV_code) %>%
  mutate(time_pop_ID=paste(time, pop_ID, treatment, sep="__")) %>% 
  ungroup() %>% 
  select(time_pop_ID, Chrom, Gene, VarFreq) %>% 
  group_by(time_pop_ID, Chrom, Gene) %>% 
  dplyr::summarize(VarFreq_max=max(VarFreq)) %>%
  merge((phaseIII_table_total_pink_top_bottom_diffGV %>% 
           separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
           filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
           filter(time==53) %>% 
           filter(GV_code %in% 
                    phaseIII_table_total_private_GVs_allo_pink_top_bottom_diffGV$GV_code) %>% 
           group_by(Gene) %>% 
           dplyr::summarize(N_GVs=length(unique(Position)), 
                            Position_min=min(Position))), all.x = T) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") 

head(phaseIII_table_total_bygene_pink_top_bottom_diffGV)

phaseIII_table_total_bygene_pink_top_bottom_diffGV<-
  cbind(phaseIII_table_total_bygene_pink_top_bottom_diffGV, 
        Gene_NumGV=factor(paste0(phaseIII_table_total_bygene_pink_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$Position_min), 
                          levels=unique(paste0(phaseIII_table_total_bygene_pink_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$Position_min)[order(phaseIII_table_total_bygene_pink_top_bottom_diffGV$Position_min)])))

phaseIII_table_total_bygene_pink_top_bottom_diffGV %>% 
  ggplot(aes(Gene_NumGV, VarFreq_max, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, shape=21, alpha=0.3, position = position_jitterdodge()) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  ggsave(paste0("plots/22_8_changeFq_Allo_phaseIII_byGene_pink.png"), width = 15, height = 8, dpi = 400)


# by gene but adding all variants in genes with at least one GV significatly different
# for instance this will add private mutations in one population in phase III

# red populations

phaseIII_table_total_bygene_red_top_bottom_Genes<-phaseIII_table_total_bygene_red_top_bottom_diffGV %>%
  ungroup() %>% 
  select(Gene) %>% unlist() %>% 
  as.vector() %>% unique()


phaseIII_table_total_bygene_red_top_bottom_diffGV<-
  phaseIII_table_total_red %>% 
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time==53) %>% 
  filter(Gene %in% phaseIII_table_total_bygene_red_top_bottom_Genes) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, 
         VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(-N_GV) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 8:29) %>% 
  ungroup() %>% 
  select(time_pop_ID, Chrom, Gene, VarFreq) %>% 
  group_by(time_pop_ID, Chrom, Gene) %>% 
  dplyr::summarize(VarFreq_max=max(VarFreq)) %>%
  merge((phaseIII_table_total_red %>% 
           filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
           filter(time==53) %>% 
           filter(Gene %in% phaseIII_table_total_bygene_red_top_bottom_Genes) %>% 
           group_by(Gene) %>% 
           dplyr::summarize(N_GVs=length(unique(Position)),Position_min=min(Position))), all.x = T) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") 


phaseIII_table_total_bygene_red_top_bottom_diffGV<-
  cbind(phaseIII_table_total_bygene_red_top_bottom_diffGV, 
        Gene_NumGV=factor(paste0(phaseIII_table_total_bygene_red_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$Position_min), 
                          levels=unique(paste0(phaseIII_table_total_bygene_red_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_red_top_bottom_diffGV$Position_min)[order(phaseIII_table_total_bygene_red_top_bottom_diffGV$Position_min)])))

phaseIII_table_total_bygene_red_top_bottom_diffGV %>% 
  ggplot(aes(Gene_NumGV, VarFreq_max, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, shape=21, alpha=0.3, position = position_jitterdodge()) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  ggsave(paste0("plots/22_8_changeFq_Allo_phaseIII_byGene_allGV_red.png"), width = 15, height = 8, dpi = 400)




# pink populations

phaseIII_table_total_bygene_pink_top_bottom_Genes<-phaseIII_table_total_bygene_pink_top_bottom_diffGV %>%
  ungroup() %>% 
  select(Gene) %>% unlist() %>% 
  as.vector() %>% unique()


phaseIII_table_total_bygene_pink_top_bottom_diffGV<-
  phaseIII_table_total_pink %>% 
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time==53) %>% 
  filter(Gene %in% phaseIII_table_total_bygene_pink_top_bottom_Genes) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, 
         VarFreq, Pvalue, Effect, Level, Gene) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID, Effect, Level, Gene) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(-N_GV) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 8:28) %>% 
  ungroup() %>% 
  select(time_pop_ID, Chrom, Gene, VarFreq) %>% 
  group_by(time_pop_ID, Chrom, Gene) %>% 
  dplyr::summarize(VarFreq_max=max(VarFreq)) %>%
  merge((phaseIII_table_total_pink %>% 
           filter(treatment %in% c("Allo_B", "Allo_T")) %>% 
           filter(time==53) %>% 
           filter(Gene %in% phaseIII_table_total_bygene_pink_top_bottom_Genes) %>% 
           group_by(Gene) %>% 
           dplyr::summarize(N_GVs=length(unique(Position)),Position_min=min(Position))), all.x = T) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") 


phaseIII_table_total_bygene_pink_top_bottom_diffGV<-
  cbind(phaseIII_table_total_bygene_pink_top_bottom_diffGV, 
        Gene_NumGV=factor(paste0(phaseIII_table_total_bygene_pink_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$Position_min), 
                          levels=unique(paste0(phaseIII_table_total_bygene_pink_top_bottom_diffGV$Gene,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$N_GVs,"_",phaseIII_table_total_bygene_pink_top_bottom_diffGV$Position_min)[order(phaseIII_table_total_bygene_pink_top_bottom_diffGV$Position_min)])))

phaseIII_table_total_bygene_pink_top_bottom_diffGV %>% 
  ggplot(aes(Gene_NumGV, VarFreq_max, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, shape=21, alpha=0.3, position = position_jitterdodge()) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  ggsave(paste0("plots/22_8_changeFq_Allo_phaseIII_byGene_allGV_pink.png"), width = 15, height = 8, dpi = 400)


#####

phaseIII_table_total_bygene_red<-
  phaseIII_table_total_red %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Effect!="intergenic_region") %>% 
  group_by(time_pop_ID, Chrom, Gene) %>% 
  dplyr::summarize(N_GVs=n(), Position_min=min(Position), VarFreq_max=max(VarFreq)) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) 

phaseIII_table_total_bygene_red<-cbind(phaseIII_table_total_bygene_red, Gene_NumGV=factor(paste0(phaseIII_table_total_bygene_red$Gene,"_",phaseIII_table_total_bygene_red$N_GVs), levels=unique(paste0(phaseIII_table_total_bygene_red$Gene,"_",phaseIII_table_total_bygene_red$N_GVs)[order(phaseIII_table_total_bygene_red$Position_min)])))


phaseIII_table_total_bygene_red %>% 
  # filter(VarFreq_max>20) %>% 
  group_by(Gene) %>% 
  dplyr::summarize(N_populations=n()) %>% 
  filter(N_populations>0) %>%
  # ungroup() %>% select(Gene) %>% unique() %>%
  # unlist() %>% as.vector() -> genes_more1pop
  ggplot(aes(N_populations)) + 
  geom_histogram(binwidth = 1) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  ggsave(paste0("plots/22_9_NumberPop_perGene_phaseIII_red.png"), width = 5, height = 4, dpi = 400)
  # ggsave(paste0("plots/22_9_NumberPop_perGene_minFq20_phaseIII_red.png"), width = 5, height = 4, dpi = 400)

# phaseIII_table_total_bygene_red %>% 
#   #filter(N_GVs==1) %>% 
#   #filter(Chrom=="I") %>% 
#   filter(Gene %in% genes_more1pop) %>% 
#   ggplot(aes(Gene_NumGV, VarFreq_max, fill=treatment)) + 
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(aes(fill = treatment), size = 2, shape=21, alpha=0.3, position = position_jitterdodge()) +
#   scale_fill_manual(values=c("red3", "cornflowerblue")) +
#   facet_grid(. ~ Chrom, scale="free", space="free") +
#   theme_classic() + 
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_line(colour = "grey70"), 
#         panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
#         axis.line = element_line(colour = "black"), 
#         #legend.position="none", 
#         axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
#         axis.text.y = element_text(size=10, colour="black"), 
#         strip.text.y = element_text(angle = 360, size=12)) +
#   ggsave(paste0("plots/22_9_changeFq_Allo_phaseIII_byGene_includesPhaseIIIGVs_red.png"), width = 15, height = 8, dpi = 400)

# pink populations:

head(phaseIII_table_total_pink)
phaseIII_table_total_bygene_pink<-
  phaseIII_table_total_pink %>% 
  filter(time==53) %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  filter(time_groupGV!="fix_parental") %>% 
  #filter(time_groupGV!="phaseIII") %>% 
  #filter(Effect!="intergenic_region") %>% 
  group_by(time_pop_ID, Chrom, Gene) %>% 
  dplyr::summarize(N_GVs=n(), Position_min=min(Position), VarFreq_max=max(VarFreq)) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(treatment %in% c("Allo_B", "Allo_T")) 

phaseIII_table_total_bygene_pink<-cbind(phaseIII_table_total_bygene_pink, Gene_NumGV=factor(paste0(phaseIII_table_total_bygene_pink$Gene,"_",phaseIII_table_total_bygene_pink$N_GVs), levels=unique(paste0(phaseIII_table_total_bygene_pink$Gene,"_",phaseIII_table_total_bygene_pink$N_GVs)[order(phaseIII_table_total_bygene_pink$Position_min)])))


phaseIII_table_total_bygene_pink %>% 
  filter(VarFreq_max>20) %>%
  group_by(Gene) %>% 
  dplyr::summarize(N_populations=n()) %>% 
  filter(N_populations>0) %>% 
  # filter(N_populations<3) %>%
  # ungroup() %>% select(Gene) %>% unique() %>%
  # unlist() %>% as.vector() -> genes_more1pop
  ggplot(aes(N_populations)) +
  geom_histogram(binwidth = 1) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  # ggsave(paste0("plots/22_9_NumberPop_perGene_phaseIII_pink.png"), width = 5, height = 4, dpi = 400)
  ggsave(paste0("plots/22_9_NumberPop_perGene_minFq20_phaseIII_pink.png"), width = 5, height = 4, dpi = 400)

phaseIII_table_total_bygene_pink %>% 
  #filter(N_GVs==1) %>% 
  #filter(Chrom=="I") %>% 
  filter(Gene %in% genes_more1pop) %>% 
  ggplot(aes(Gene_NumGV, VarFreq_max, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, shape=21, alpha=0.3, position = position_jitterdodge()) +
  scale_fill_manual(values=c("red3", "cornflowerblue")) +
  facet_grid(. ~ Chrom, scale="free", space="free") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(size=10, colour="black"), 
        strip.text.y = element_text(angle = 360, size=12)) +
  ggsave(paste0("plots/22_9_changeFq_Allo_phaseIII_byGene_includesPhaseIIIGVs_pink.png"), width = 18, height = 5, dpi = 400)





# PCA
library(ggfortify)
library(factoextra)
  
# red populations:
phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
             (time!=(-16) | VarFreq>70) &
             (time!=(-10) | VarFreq>70) &
             time %in% c(-16,-10,-5,0,53) &
             parental == "P3R" ) 
  
phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))
  
  
# allopatric populations:
phaseIII_table_total_PCA_red<-  phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:27) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#D95F02")
# colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
  
PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))

  
autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=F, 
         label.size = 5, 
         frame=T)  +
  geom_point(size=4, aes(colour=treatment)) +
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.title = element_text(size=24, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) #+
  ggsave("plots/23_PCA_total_allo_red_vector.svg", width = 7, height = 5)
# ggsave("plots/23_PCA_total_allo_red.png", width = 7, height = 5, dpi = 400)
# ggsave("plots/23_PCA_ancestral_allo_red.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() #+
  # ggsave("plots/23_PCA_total_allo_hist_red.png", width = 6, height = 6, dpi = 400)
# ggsave("plots/23_PCA_ancestral_allo_hist_red.png", width = 6, height = 6, dpi = 400)

  
# testing significants with GLMs:

library("lme4")
  
res.ind <- get_pca_ind(PCA_data)

mod_allo<-lm(Dim.1 ~ treatment, 
             data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_red[, 1:3]) %>% 
               filter(time==53))
res_mod<-summary(mod_allo)
res_mod

mod_allo<-lm(Dim.2 ~ treatment, 
             data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_red[, 1:3]) %>% 
               filter(time==53))
res_mod<-summary(mod_allo)
res_mod



  
  
# test for significance of eigenvalues:

table_test_pvalue<-phaseIII_table_total_PCA_red %>% 
  select(-treatment, -time, -pop_ID) 

Y0<-as.matrix(table_test_pvalue)
rownames(Y0)<-(phaseIII_table_total_PCA_red %>% 
                                 mutate(populationID=paste(treatment, time, pop_ID, sep="__")) %>% 
                                 select(populationID) %>% 
                                 unlist() %>% 
                                 as.vector())
alpha <- 0.05
N.Boot <- 1000
n <- nrow(Y0) ; n   
p <- ncol(Y0) ; p
M <- min(n - 1, p) ; M
# Variables are centred, but not scaled
Y <- scale(Y0, center = TRUE, scale = FALSE)
# Singular value decomposition
svd.Y <- svd(Y)
lambda <- svd.Y$d
T <- {lambda^2/rev(cumsum(rev(lambda^2)))} [1:(M - 1)]
# Sampling
T.boot <- matrix(NA, nrow = N.Boot, ncol = M - 1)
for(Component in 1:(M - 1)){
  print(Component)
  K <- Component - 1
  for(Boot in 1:N.Boot){
    E.boot <- matrix(rnorm((n - 1 - K)*(p - K)), nrow = n - 1 - K, ncol = p - K)
    lambda.boot <- svd(E.boot)$d
    T.boot[Boot, Component] <- lambda.boot[1]^2 / sum(lambda.boot^2)
  }
  if (mean(T.boot[, Component] > T[Component]) > alpha){break}
}
# Results
N.Boot
lambda^2
T
colMeans(T.boot > matrix(rep(T, N.Boot), nrow = N.Boot, byrow = TRUE)) # p-values

library(AssocTests)
phaseIII_table_total_PCA_red %>% 
  head()

table_test_pvalue<-phaseIII_table_total_PCA_red %>% 
  select(-treatment, -time, -pop_ID) 

rownames(table_test_pvalue/100)

table_test_pvalue2<-Y0<-as.matrix(table_test_pvalue/100)
table_test_pvalue2[,1]

r<-dim(table_test_pvalue2)[1]
c<-dim(table_test_pvalue2)[2]
total_samples<-10

m00 <- matrix(0, r, c)
for(i in 1:r){
  for(j in 1:c){
    m00[i, j] <- rbinom(1, size = 1, prob=table_test_pvalue2[i,j])
  }
}

for (sample_used in seq(2,total_samples)){
  print(sample_used)
  m01 <- matrix(0, r, c)
  for(i in 1:r){
    for(j in 1:c){
      m01[i, j] <- rbinom(1, size = 1, prob=table_test_pvalue2[i,j])
    }
  }
  m00<-rbind(m00, m01)
}


head(m00)
dim(m00)
eigenstratG.eg <- m00
getwd()
write.table(eigenstratG.eg, file = "eigenstratG.eg.txt", quote = FALSE,
            sep = "", row.names = FALSE, col.names = FALSE)
a<-eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
           outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
           rm.subject.index = NULL, miss.val = 9, num.splits = 10,
           topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
           iter.outlier = 5, sigma.thresh = 6)

head(a)
str(a)
plot(a$eigenvectors[,1], a$eigenvectors[,2])
a$TW.stat

eigenstrat<-function(geno,maxMis=0,minMaf=0.01){                 
  ## geno: snp x ind matrix of genotypes \in 0,1,2
  ##maxMis maximum allowed missing genotypes for a site
  
  nMis <- rowSums(is.na(geno))
  freq <- rowMeans(geno,na.rm=T)/2               # get allele frequency 
  keep <- freq>minMaf&freq<(1-minMaf) & nMis<=maxMis         # remove sites with non-polymorphic data
  freq<-freq[keep]
  geno<-geno[keep,]
  snp<-nrow(geno)                           #number of snps used in analysis
  ind<-ncol(geno)                           #number of individuals used in analuysis
  M <- (geno-freq*2)/sqrt(freq*(1-freq))       #normalize the genotype matrix
  M[is.na(M)]<-0
  X<-t(M)%*%M                               #get the (almost) covariance matrix
  X<-X/(sum(diag(X))/(snp-1))
  E<-eigen(X)
  
  mu<-(sqrt(snp-1)+sqrt(ind))^2/snp         #for testing significance (assuming no LD!)
  sigma<-(sqrt(snp-1)+sqrt(ind))/snp*(1/sqrt(snp-1)+1/sqrt(ind))^(1/3)
  E$TW<-(E$values[1]*ind/sum(E$values)-mu)/sigma
  E$mu<-mu
  E$sigma<-sigma
  E$nSNP <- nrow(geno)
  E$nInd <- ncol(geno)
  class(E)<-"eigenstrat"
  E
}



plot.eigenstrat<-function(x,col=1,...)
  plot(x$vectors[,1:2],col=col,...)

print.eigenstrat<-function(x)
  cat("statistic",x$TW,"\n")

ind<-c(100,100,100)
snp<-10000
freq<-cbind(runif(snp),runif(snp),runif(snp))
l<-lapply(1:length(ind),function(x) matrix(rbinom(snp*ind[x],2,freq[,x]),snp))
geno<-do.call(cbind,l)
e<-eigenstrat(geno)
str(e)
e$sigma
plot(e,col=rep(1:length(ind),ind),xlab="PC1",ylab="PC2")

plot.eigenstrat(e)
print.eigenstrat(e)
e

# names(PCA_data)
# dim(PCA_data$rotation)
# head(PCA_data$rotation)
# PCA_data$rotation[695,]
# head(PCA_data$center)
# str(PCA_data$center)
# PCA_data$center[695]
# hist(PCA_data$center[PCA_data$center>20])
# PCA_data$sdev
# head(PCA_data$x)
# dim(PCA_data$x)
# plot(PCA_data$x[,1],PCA_data$x[,2])
# which(row.names(PCA_data$rotation)=="II_1718756_A")

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = F ) #+ 
  ggsave("plots/23_PCA_total_allo_GVs2_red.png", width = 10, height = 10, dpi = 400)
# ggsave("plots/23_PCA_ancestral_allo_GVs2_red.png", width = 10, height = 10, dpi = 400)

  
# Other treatments:
phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))
  
autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=F, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  ggsave("plots/23_PCA_total_red.png", width = 7, height = 5, dpi = 400)
  # ggsave("plots/23_PCA_ancestral_red.png", width = 7, height = 5, dpi = 400)

  
fviz_eig(PCA_data) +
    theme_classic() +
    # ggsave("plots/23_PCA_total_hist_red.png", width = 6, height = 6, dpi = 400)
    ggsave("plots/23_PCA_ancestral_hist_red.png", width = 6, height = 6, dpi = 400)
  
fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_total_GVs_red.png", width = 10, height = 10, dpi = 400)
ggsave("plots/23_PCA_ancestral_GVs2_red.png", width = 10, height = 10, dpi = 400)


# Significant test using LM:

library("lme4")

res.ind <- get_pca_ind(PCA_data)
res.ind$coord %>%
  cbind(phaseIII_table_total_PCA_red[, 1:3]) %>% 
  filter(time==53) %>% 
  ggplot(aes(Dim.1, Dim.2, colour=treatment)) +
  geom_point()+
  facet_grid(treatment ~ .)

mod_allpop<-lm(Dim.1 ~ treatment, 
             data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_red[, 1:3]) %>% 
               filter(time==53))
res_mod<-summary(mod_allpop)
res_mod

mod_allpop<-lm(Dim.2 ~ treatment, 
             data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_red[, 1:3]) %>% 
               filter(time==53))
res_mod<-summary(mod_allpop)
res_mod

mod_allpop <- lm(cbind(Dim.1, Dim.2) ~ treatment, 
           data = res.ind$coord %>% 
             cbind(phaseIII_table_total_PCA_red[, 1:3]) %>% 
             filter(time==53))
summary(mod_allpop)

library(car)
Anova(mod_allpop)


# test for significance of igenvalues:

table_test_pvalue<-phaseIII_table_total_PCA_red %>% 
  select(-treatment, -time, -pop_ID) 

head(phaseIII_table_total_PCA_red)
phaseIII_table_total_PCA_red %>% 
  mutate(populationID=paste(treatment, time, pop_ID, sep="__")) %>% 
  select(populationID) %>% 
  unlist() %>% 
  as.vector() 

Y0<-as.matrix(table_test_pvalue)
rownames(Y0)<-(phaseIII_table_total_PCA_red %>% 
                                 mutate(populationID=paste(treatment, time, pop_ID, sep="__")) %>% 
                                 select(populationID) %>% 
                                 unlist() %>% 
                                 as.vector())
alpha <- 0.05
N.Boot <- 1000
n <- nrow(Y0) ; n   
p <- ncol(Y0) ; p
M <- min(n - 1, p) ; M
# Variables are centred, but not scaled
Y <- scale(Y0, center = TRUE, scale = FALSE)
# Singular value decomposition
svd.Y <- svd(Y)
lambda <- svd.Y$d
T <- {lambda^2/rev(cumsum(rev(lambda^2)))} [1:(M - 1)]
# Sampling
T.boot <- matrix(NA, nrow = N.Boot, ncol = M - 1)
for(Component in 1:(M - 1)){
  print(Component)
  K <- Component - 1
  for(Boot in 1:N.Boot){
    E.boot <- matrix(rnorm((n - 1 - K)*(p - K)), nrow = n - 1 - K, ncol = p - K)
    lambda.boot <- svd(E.boot)$d
    T.boot[Boot, Component] <- lambda.boot[1]^2 / sum(lambda.boot^2)
  }
  if (mean(T.boot[, Component] > T[Component]) > alpha){break}
}
# Results
N.Boot
lambda^2
T
colMeans(T.boot > matrix(rep(T, N.Boot), nrow = N.Boot, byrow = TRUE)) # p-values

phaseIII_table_total_PCA_red %>% 
  head()




# pink populations:

library(ggfortify)
library(factoextra)

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))

# allopatric populations:
phaseIII_table_total_PCA_pink<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>%
  # filter(time_groupGV!="phaseIII") %>%
  filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>% 
  # filter(time %in% c(53)) %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:26) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3P")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#D95F02")
# colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_pink, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  # ggsave("plots/23_PCA_total_allo_pink.png", width = 7, height = 5, dpi = 400)
  ggsave("plots/23_PCA_ancestral_allo_pink.png", width = 7, height = 5, dpi = 400)

fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_total_allo_hist_pink.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_allo_hist_pink.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_total_allo_GVs2_pink.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_allo_GVs2_pink.png", width = 10, height = 10, dpi = 400)


# testing significants with GLMs:

library("lme4")

res.ind <- get_pca_ind(PCA_data)

mod_allo<-lm(Dim.1 ~ treatment, 
             data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_pink[, 1:3]) %>% 
               filter(time==53))
res_mod<-summary(mod_allo)
res_mod

mod_allo<-lm(Dim.2 ~ treatment, 
             data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_pink[, 1:3]) %>% 
               filter(time==53))
res_mod<-summary(mod_allo)
res_mod




# Other treatments:

phaseIII_table_total_PCA_pink<-phaseIII_table_total_pink %>%
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  #filter(time %in% c(53)) %>%
  #filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:70) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3P", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_pink, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  # ggsave("plots/23_PCA_total_pink.png", width = 7, height = 5, dpi = 400)
  ggsave("plots/23_PCA_ancestral_pink.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_total_hist_pink.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_hist_pink.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_total_GVs_pink.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_GVs2_pink.png", width = 10, height = 10, dpi = 400)


# Significant test using LM:

library("lme4")

res.ind <- get_pca_ind(PCA_data)
res.ind$coord %>%
  cbind(phaseIII_table_total_PCA_pink[, 1:3]) %>% 
  filter(time==53) %>% 
  ggplot(aes(Dim.1, Dim.2, colour=treatment)) +
  geom_point()+
  facet_grid(treatment ~ .)

mod_allpop<-lm(Dim.1 ~ treatment, 
               data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_pink[, 1:3]) %>% 
                 filter(time==53))
res_mod<-summary(mod_allpop)
res_mod

mod_allpop<-lm(Dim.2 ~ treatment, 
               data=res.ind$coord %>% cbind(phaseIII_table_total_PCA_pink[, 1:3]) %>% 
                 filter(time==53))
res_mod<-summary(mod_allpop)
res_mod

mod_allpop <- lm(cbind(Dim.1, Dim.2) ~ treatment, 
                 data = res.ind$coord %>% 
                   cbind(phaseIII_table_total_PCA_pink[, 1:3]) %>% 
                   filter(time==53))
summary(mod_allpop)
library(car)
Anova(mod_allpop)


# Allo Vs other treatments independently:

# red populations:
# Parapatric:
phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>%
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>%
  # filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  # gather("time_pop_ID", "VarFreq", 5:27) %>% 
  gather("time_pop_ID", "VarFreq", 5:49) %>%
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#D95F02", "#1B9E77", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  geom_point(size=4, aes(colour=treatment), alpha=0.8) +
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.title = element_text(size=24, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  ggsave("plots/23_PCA_para_allo_red_vector.svg", width = 7, height = 5)
  # ggsave("plots/23_PCA_para_allo_red.png", width = 7, height = 5, dpi = 400)
  # ggsave("plots/23_PCA_ancestral_para_allo_red.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_para_allo_hist_red.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_para_allo_hist_red.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_para_allo_GVs_red.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_para_allo_GVs_red.png", width = 10, height = 10, dpi = 400)


# Only Parapatric:
phaseIII_table_total_PCA_red<-phaseIII_table_total_red %>%
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:27) %>% 
  # gather("time_pop_ID", "VarFreq", 5:49) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("#D95F02", "#1B9E77", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))

# row.names(phaseIII_table_total_PCA_red)<-phaseIII_table_total_PCA_red$pop_ID
row.names(phaseIII_table_total_PCA_red)<-phaseIII_table_total_PCA_red %>% 
  separate(pop_ID, c("pop", "treat", "frac")) %>% 
  mutate(pop_ID=paste(pop,frac,sep="_")) %>% 
  ungroup() %>% 
  select(pop_ID) %>% 
  unlist() %>% 
  as.vector()

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = T , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=F, label.size = 5, 
         frame=F) + 
  geom_point(size=4, aes(colour=treatment), alpha=0.8) +
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.title = element_text(size=24, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) #+
  # ggsave("plots/23_PCA_para_red_vector.svg", width = 7, height = 5)

fviz_eig(PCA_data) +
  theme_classic() 

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_para_allo_GVs_red.png", width = 10, height = 10, dpi = 400)
  # ggsave("plots/23_PCA_para_GVs2_red.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_para_GVs2_red.svg", width = 10, height = 10)





#####

# Allo Vs other treatments independently:
# Sympatry:
phaseIII_table_total_PCA_red<-
  phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM")) %>% 
  # filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:37) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#D95F02", "#E6AB02")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))



autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  geom_point(size=4, aes(colour=treatment), alpha=0.8) +
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.title = element_text(size=24, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  ggsave("plots/23_PCA_sym_allo_red_vector.svg", width = 7, height = 5)
# ggsave("plots/23_PCA_sym_allo_red.png", width = 7, height = 5, dpi = 400)
# ggsave("plots/23_PCA_ancestral_sym_allo_red.png", width = 7, height = 5, dpi = 400)



autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) #+
  # ggsave("plots/23_PCA_sym_allo_red.png", width = 7, height = 5, dpi = 400)
  ggsave("plots/23_PCA_ancestral_sym_allo_red.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_sym_allo_hist_red.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_sym_allo_hist_red.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_sym_allo_GVs_red.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_sym_allo_GVs2_red.png", width = 10, height = 10, dpi = 400)

#####




# Allo Vs other treatments independently:
# Local Mating:
phaseIII_table_total_PCA_red<-
  phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "LM")) %>% 
  # filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:38) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#66A61E", "#D95F02")

PCA_data<-prcomp(phaseIII_table_total_PCA_red %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_red, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  geom_point(size=4, aes(colour=treatment), alpha=0.8) +
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.title = element_text(size=24, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  ggsave("plots/23_PCA_lm_allo_red_vector.svg", width = 7, height = 5)
  # ggsave("plots/23_PCA_lm_allo_red.png", width = 7, height = 5, dpi = 400)
  # ggsave("plots/23_PCA_ancestral_lm_allo_red.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  ggsave("plots/23_PCA_lm_allo_hist_red.png", width = 6, height = 6, dpi = 400)
# ggsave("plots/23_PCA_ancestral_lm_allo_hist_red.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F ,
             repel = F ) + 
  ggsave("plots/23_PCA_lm_allo_GVs_red.png", width = 10, height = 10, dpi = 400)
# ggsave("plots/23_PCA_ancestral_lm_allo_GVs2_red.png", width = 10, height = 10, dpi = 400)


####


# pink populations:


# Parapatric:
phaseIII_table_total_PCA_pink<-
  phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:47) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_pink, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  # ggsave("plots/23_PCA_para_allo_pink.png", width = 7, height = 5, dpi = 400)
  ggsave("plots/23_PCA_ancestral_para_allo_pink.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_para_allo_hist_pink.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_para_allo_hist_pink.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_para_allo_GVs2_pink.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_para_allo_GVs_pink.png", width = 10, height = 10, dpi = 400)



#####

# Allo Vs other treatments independently:
# Sympatry:
phaseIII_table_total_PCA_pink<-
  phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM")) %>% 
  filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:36) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_pink, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  # ggsave("plots/23_PCA_sym_allo_pink.png", width = 7, height = 5, dpi = 400)
  ggsave("plots/23_PCA_ancestral_sym_allo_pink.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_sym_allo_hist_pink.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_sym_allo_hist_pink.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_sym_allo_GVs2_pink.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_sym_allo_GVs_pink.png", width = 10, height = 10, dpi = 400)


#####


# Allo Vs other treatments independently:
# Local Mating:
phaseIII_table_total_PCA_pink<-
  phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "LM")) %>% 
  filter(time_groupGV!="phaseIII") %>%
  filter(time %in% c(0,53)) %>%
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:36) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  group_by(treatment, time, pop_ID) %>% 
  select(GV_code, treatment, time, pop_ID, VarFreq) %>% 
  spread(GV_code, VarFreq) %>% 
  ungroup() 

colours_PCA<-c("red3", "cornflowerblue", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")

PCA_data<-prcomp(phaseIII_table_total_PCA_pink %>% select(-treatment, -time, -pop_ID))

autoplot(PCA_data, 
         data=phaseIII_table_total_PCA_pink, 
         label = F , 
         colour = 'treatment', 
         fill =  'treatment', 
         shape=T, label.size = 5, 
         frame=T) + 
  scale_colour_manual(values=colours_PCA) + 
  scale_fill_manual(values=colours_PCA) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #legend.position="none", 
        axis.text.x=element_text(size=13, colour="black"), 
        axis.text.y = element_text(colour="black", size=13)) +
  # ggsave("plots/23_PCA_lm_allo_pink.png", width = 7, height = 5, dpi = 400)
  ggsave("plots/23_PCA_ancestral_lm_allo_pink.png", width = 7, height = 5, dpi = 400)


fviz_eig(PCA_data) +
  theme_classic() +
  # ggsave("plots/23_PCA_lm_allo_hist_pink.png", width = 6, height = 6, dpi = 400)
  ggsave("plots/23_PCA_ancestral_lm_allo_hist_pink.png", width = 6, height = 6, dpi = 400)

fviz_pca_var(PCA_data,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = F ,
             repel = F ) + 
  # ggsave("plots/23_PCA_lm_allo_GVs2_pink.png", width = 10, height = 10, dpi = 400)
  ggsave("plots/23_PCA_ancestral_lm_allo_GVs2_pink.png", width = 10, height = 10, dpi = 400)



#####

# comparing between top and bottom differentiated GVs in all treatments

phaseIII_table_total_red %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_pop_ID, time_groupGV,
         VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_pop_ID, time_groupGV) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(-N_GV) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  group_by(time_pop_ID) %>% 
  filter(GV_code %in% c("II_1741521_T", "II_1718756_A")) %>% 
  select(time_pop_ID, GV_code, VarFreq) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(time %in% c(53)) %>% 
  spread(GV_code, VarFreq, fill=(-100)) %>% 
  mutate(treatment=factor(treatment, levels=c("Allo_T", "Allo_B", "Para_T", "Para_B", "LM", "SYM"))) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  # filter(treatment %in% c("Para_T", "Para_B")) %>%
  filter(treatment %in% c("LM", "SYM")) %>%
  ggplot(aes(II_1718756_A/100, II_1741521_T/100)) + #, colour=treatment)) +
  geom_abline(intercept = 1, slope = -1, alpha=0.6) +
  geom_jitter(size=2, alpha=0.7, width = 0.015, height = 0.015) + #geom_point(size=2, alpha=0.4) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  facet_grid(. ~ treatment) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey85"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom",
        strip.background = element_blank(), 
        # legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
# ggsave(paste0("plots/22_7_II_1741521_T_II_1718756_A_phaseIII_red.png"), width = 5, height = 6, dpi = 400)
# ggsave(paste0("plots/22_7_II_1741521_T_II_1718756_A_2_phaseIII_red.png"), width = 4, height = 10, dpi = 400)
ggsave(paste0("plots/22_7_II_1741521_T_II_1718756_A_2_phaseIII_red_3.png"), width = 5, height = 3, dpi = 400)


phaseIII_table_total_pink %>% 
  mutate(time_pop_ID=paste0(time,"__",pop_treatment,"__",treatment)) %>% 
  select(Chrom, Position, GV_code, time_pop_ID, time_groupGV,
         VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_pop_ID, time_groupGV) %>% 
  dplyr::summarize(N_GV=n(), VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(-N_GV) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:85) %>% 
  group_by(time_pop_ID) %>% 
  filter(GV_code %in% c("I_2319886_T", "I_2319922_T")) %>% 
  select(time_pop_ID, GV_code, VarFreq) %>% 
  separate(time_pop_ID, c("time", "pop_ID", "treatment"), "__") %>% 
  filter(time %in% c(0,53)) %>% 
  spread(GV_code, VarFreq, fill=(-100)) %>% 
  mutate(treatment=factor(treatment, levels=c("Allo_T", "Allo_B", "Para_T", "Para_B", "LM", "SYM"))) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B")) %>%
  # filter(treatment %in% c("Para_T", "Para_B")) %>%
  # filter(treatment %in% c("LM", "SYM")) %>%
  filter(is.na(treatment)) %>%
  ggplot(aes(I_2319886_T/100, I_2319922_T/100)) + #, colour=treatment)) +
  geom_abline(intercept = 1, slope = -1, alpha=0.6) +
  geom_jitter(size=2, alpha=0.7, width = 0.015, height = 0.015) + #geom_point(size=2, alpha=0.4) +
  ylim(c(-0.02,1)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  facet_grid(. ~ treatment) +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey85"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom",
        strip.background = element_blank(), 
        # legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("plots/22_7_I_2319886_T_I_2319922_T_phaseIII_pink.png"), width = 5, height = 6, dpi = 400)
  ggsave(paste0("plots/22_7_I_2319886_T_I_2319922_T_phaseIII_pink_1.png"), width = 5, height = 3, dpi = 400)
  




#####

# testing for linkage (LD) beyween pairs of GVS:

# test example: just to figure out how the model works:
N<-150000
s=0.1
x1<-1/(2*N)
x2<-0
x3<-(0.5)-1/(2*N)
x4<-0.5
r<-0.001
dis<-c()
h<-c()
p2_final<-c()

for (r in seq(0.0001,0.4,0.0001)) {
  dis<-c(dis, (r/s))
  while ((x1+x2) < (1-(1/(2*N))) ){
    # print(c(x1, x2))
    # print(x1 + x2)
    p1<- (x1 + x2)
    q1<- (1-p1)
    wbar1 <- (1-(q1*s*0.5))
    wbar3 <- (1-(p1*s*0.5)-(q1*s))
    wbar <- (1-(q1*s))
    rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
    x1 <- ((x1*wbar1)-rWD)/wbar
    x2 <- ((x2*wbar1)+rWD)/wbar
    x3 <- ((x3*wbar3)+rWD)/wbar
    x4 <- ((x4*wbar3)-rWD)/wbar
  }
  p2 <- x1 + x3
  q2 <- 1- p2
  h<-c(h, 2*p2*q2/0.5)
  p2_final <- c(p2_final, p2)
  x1<-1/(2*N)
  x2<-0
  x3<-(0.5)-1/(2*N)
  x4<-0.5
}

# table<-data.frame(dis, h, p2_final, sel=s, Ne=N)
table<-rbind(table, data.frame(dis, h, p2_final, sel=s, Ne=N))

table %>% 
  ggplot(aes(dis,p2_final, colour=as.factor(Ne))) +
  geom_line() +
  xlim(c(0,1)) +
  facet_grid(as.factor(sel) ~ .)


# testing red pair:
# first botton GV is under selection and top is neutral:
# Selection in a bottom population:

# considering initial parental frequencies for the two GVs, I produce a list or posible haplotype genotypes frequencies:
ini_x1v<-c()
ini_x2v<-c()
ini_x3v<-c()
ini_x4v<-c()
for(A1 in seq(0.02,1, 0.02)) {
  ini_x1<-A1
  ini_x2<-0.5-ini_x1
  ini_x3<-0.8-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1v<-c(ini_x1v, ini_x1)
  ini_x2v<-c(ini_x2v, ini_x2)
  ini_x3v<-c(ini_x3v, ini_x3)
  ini_x4v<-c(ini_x4v, ini_x4)
}
A1v<-data.frame(ini_x1v, ini_x2v, ini_x3v, ini_x4v) %>% 
  filter(ini_x1v>0 & ini_x2v>0 & ini_x3v>0 & ini_x4v>0) %>% 
  select(ini_x1v) %>%
  arrange(ini_x1v) %>% 
  unlist () %>% 
  as.vector() 
# ggplot(aes(ini_x1v, ini_x2v, colour="red")) +
#   geom_line() +
#   geom_line(aes(ini_x1v, ini_x3v, colour="blue")) +
#   geom_line(aes(ini_x1v, ini_x4v, colour="black")) 

N<-15000
dis<-c()
p2_final<-c()
A1_value<-c()
pop<-c()
s=0.05
# r<-0.0001
# A1<-0.25

for(A1 in A1v) {
  ini_x1<-A1
  ini_x2<-0.5-ini_x1
  ini_x3<-0.8-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1T<-ini_x4
  ini_x2T<-ini_x3
  ini_x3T<-ini_x2
  ini_x4T<-ini_x1
  # Selection in a bottom population:
  x1<-ini_x1; x2<-ini_x2; x3<-ini_x3; x4<-ini_x4
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Bottom pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, p2)
    x1<-ini_x1
    x2<-ini_x2
    x3<-ini_x3
    x4<-ini_x4
  }
  # Selection in a top population:
  x1<-ini_x1T; x2<-ini_x2T; x3<-ini_x3T; x4<-ini_x4T
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Top pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, q2)
    x1<-ini_x1T
    x2<-ini_x2T
    x3<-ini_x3T
    x4<-ini_x4T
  }
}

# table_bottomSelectedlocus<-data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop)
table_bottomSelectedlocus<-rbind(table_bottomSelectedlocus, data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop))

table_bottomSelectedlocus %>% 
  ggplot(aes(dis,p2_final, colour=A1, group=A1)) +
  geom_line() +
  ylim(c(0,1)) +
  # xlim(c(0,0.1)) +
  facet_grid(as.factor(sel) ~ population) + 
  labs(y = "Final frequency neutral locus", 
       x = "Distance to selected locus (r)", 
       colour= "Initial Fq. A1B1 Genotype") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom",
        # legend.position="none", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
ggsave(paste0("plots/22_7_LDsimulation_BottomSel_II_V1_1741521_V2_1718756_phIII_red.png"), width = 7, height = 6, dpi = 400)


# Now top GV is under selection and bottom is neutral:
ini_x1v<-c()
ini_x2v<-c()
ini_x3v<-c()
ini_x4v<-c()
for(A1 in seq(0.02,1, 0.02)) {
  ini_x1<-A1
  ini_x2<-0.2-ini_x1
  ini_x3<-0.5-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1v<-c(ini_x1v, ini_x1)
  ini_x2v<-c(ini_x2v, ini_x2)
  ini_x3v<-c(ini_x3v, ini_x3)
  ini_x4v<-c(ini_x4v, ini_x4)
}


A1v<-
  data.frame(ini_x1v, ini_x2v, ini_x3v, ini_x4v) %>% 
  filter(ini_x1v>0 & ini_x2v>0 & ini_x3v>0 & ini_x4v>0) %>%
  select(ini_x1v) %>%
  arrange(ini_x1v) %>%
  unlist () %>%
  as.vector()
# ggplot(aes(ini_x1v, ini_x2v, colour="red")) +
#   geom_line() +
#   geom_line(aes(ini_x1v, ini_x3v, colour="blue")) +
#   geom_line(aes(ini_x1v, ini_x4v, colour="black")) 

N<-15000
s=0.05
dis<-c()
p2_final<-c()
A1_value<-c()
pop<-c()
# r<-0.0001
# A1<-0.25

for(A1 in A1v) {
  ini_x1<-A1
  ini_x2<-0.2-ini_x1
  ini_x3<-0.5-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1T<-ini_x4
  ini_x2T<-ini_x3
  ini_x3T<-ini_x2
  ini_x4T<-ini_x1
  # Selection in a top population:
  x1<-ini_x1; x2<-ini_x2; x3<-ini_x3; x4<-ini_x4
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Top pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, p2)
    x1<-ini_x1
    x2<-ini_x2
    x3<-ini_x3
    x4<-ini_x4
  }
  # Selection in a bottom population:
  x1<-ini_x1T; x2<-ini_x2T; x3<-ini_x3T; x4<-ini_x4T
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Bottom pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, q2)
    x1<-ini_x1T
    x2<-ini_x2T
    x3<-ini_x3T
    x4<-ini_x4T
  }
}

# table_topSelectedlocus<-data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop)
table_topSelectedlocus<-rbind(table_topSelectedlocus, data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop))


table_topSelectedlocus %>% 
  mutate(population=factor(population, levels=c("Top pop.", "Bottom pop."))) %>% 
  ggplot(aes(dis,p2_final, colour=A1, group=A1)) +
  geom_line() +
  ylim(c(0,1)) +
  # xlim(c(0,0.1)) +
  facet_grid(as.factor(sel) ~ population) + 
  labs(y = "Final frequency neutral locus", 
       x = "Distance to selected locus (r)", 
       colour= "Initial Fq. A1B1 Genotype") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom",
        # legend.position="none", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave(paste0("plots/22_7_LDsimulation_TopSel_II_V1_1741521_V2_1718756_phIII_red.png"), width = 7, height = 6, dpi = 400)

###



# testing pink pair:
# first botton GV is under selection and top is neutral:
# Selection in a bottom population:

# considering initial parental frequencies for the two GVs, I produce a list or posible haplotype genotypes frequencies:
ini_x1v<-c()
ini_x2v<-c()
ini_x3v<-c()
ini_x4v<-c()
for(A1 in seq(0.01,1, 0.01)) {
  ini_x1<-A1
  ini_x2<-0.6-ini_x1
  ini_x3<-0.9-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1v<-c(ini_x1v, ini_x1)
  ini_x2v<-c(ini_x2v, ini_x2)
  ini_x3v<-c(ini_x3v, ini_x3)
  ini_x4v<-c(ini_x4v, ini_x4)
}
A1v<-
  data.frame(ini_x1v, ini_x2v, ini_x3v, ini_x4v) %>% 
  filter(ini_x1v>0 & ini_x2v>0 & ini_x3v>0 & ini_x4v>0) %>%
  select(ini_x1v) %>%
  arrange(ini_x1v) %>%
  unlist () %>%
  as.vector()
# ggplot(aes(ini_x1v, ini_x2v, colour="pink")) +
#   geom_line() +
#   geom_line(aes(ini_x1v, ini_x3v, colour="blue")) +
#   geom_line(aes(ini_x1v, ini_x4v, colour="black"))

N<-15000
dis<-c()
p2_final<-c()
A1_value<-c()
pop<-c()
s=0.05
# r<-0.0001
# A1<-0.25

for(A1 in A1v) {
  ini_x1<-A1
  ini_x2<-0.6-ini_x1
  ini_x3<-0.9-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1T<-ini_x4
  ini_x2T<-ini_x3
  ini_x3T<-ini_x2
  ini_x4T<-ini_x1
  # Selection in a bottom population:
  x1<-ini_x1; x2<-ini_x2; x3<-ini_x3; x4<-ini_x4
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Bottom pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, p2)
    x1<-ini_x1
    x2<-ini_x2
    x3<-ini_x3
    x4<-ini_x4
  }
  # Selection in a top population:
  x1<-ini_x1T; x2<-ini_x2T; x3<-ini_x3T; x4<-ini_x4T
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Top pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, q2)
    x1<-ini_x1T
    x2<-ini_x2T
    x3<-ini_x3T
    x4<-ini_x4T
  }
}

# table_bottomSelectedlocus<-data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop)
table_bottomSelectedlocus<-rbind(table_bottomSelectedlocus, data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop))

table_bottomSelectedlocus %>% 
  ggplot(aes(dis,p2_final, colour=A1, group=A1)) +
  geom_line() +
  ylim(c(0,1)) +
  # xlim(c(0,0.1)) +
  facet_grid(as.factor(sel) ~ population) + 
  labs(y = "Final frequency neutral locus", 
       x = "Distance to selected locus (r)", 
       colour= "Initial Fq. A1B1 Genotype") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom",
        # legend.position="none", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  ggsave(paste0("plots/22_7_LDsimulation_BottomSel_I_V1_2319886_V2_2319922_phIII_pink.png"), width = 7, height = 6, dpi = 400)


# Now top GV is under selection and bottom is neutral:
ini_x1v<-c()
ini_x2v<-c()
ini_x3v<-c()
ini_x4v<-c()
for(A1 in seq(0.01,1, 0.01)) {
  ini_x1<-A1
  ini_x2<-0.1-ini_x1
  ini_x3<-0.4-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1v<-c(ini_x1v, ini_x1)
  ini_x2v<-c(ini_x2v, ini_x2)
  ini_x3v<-c(ini_x3v, ini_x3)
  ini_x4v<-c(ini_x4v, ini_x4)
}


A1v<-
  data.frame(ini_x1v, ini_x2v, ini_x3v, ini_x4v) %>% 
  filter(ini_x1v>0 & ini_x2v>0 & ini_x3v>0 & ini_x4v>0) %>%
  select(ini_x1v) %>%
  arrange(ini_x1v) %>%
  unlist () %>%
  as.vector()
# ggplot(aes(ini_x1v, ini_x2v, colour="pink")) +
#   geom_line() +
#   geom_line(aes(ini_x1v, ini_x3v, colour="blue")) +
#   geom_line(aes(ini_x1v, ini_x4v, colour="black"))

N<-15000
s=0.05
dis<-c()
p2_final<-c()
A1_value<-c()
pop<-c()
# r<-0.0001
# A1<-0.25

for(A1 in A1v) {
  ini_x1<-A1
  ini_x2<-0.1-ini_x1
  ini_x3<-0.4-ini_x1
  ini_x4<-1-ini_x1-ini_x2-ini_x3
  ini_x1T<-ini_x4
  ini_x2T<-ini_x3
  ini_x3T<-ini_x2
  ini_x4T<-ini_x1
  # Selection in a top population:
  x1<-ini_x1; x2<-ini_x2; x3<-ini_x3; x4<-ini_x4
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Top pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, p2)
    x1<-ini_x1
    x2<-ini_x2
    x3<-ini_x3
    x4<-ini_x4
  }
  # Selection in a bottom population:
  x1<-ini_x1T; x2<-ini_x2T; x3<-ini_x3T; x4<-ini_x4T
  for (r in seq(0.001,0.5,0.001)) {
    pop<-c(pop,"Bottom pop.")
    dis<-c(dis, r)
    A1_value<-c(A1_value, A1)
    while ((x1+x2) < 0.9999){
      # for (n in seq(1,153)) {
      # print(n)
      # print(c(x1, x2))
      # print(x1 + x2)
      p1<- (x1 + x2)
      q1<- (1-p1)
      wbar1 <- (1-(q1*s*0.5))
      wbar3 <- (1-(p1*s*0.5)-(q1*s))
      wbar <- (1-(q1*s))
      rWD <- (r*((1-s)/2)*((x1*x4)-(x2*x3)))
      x1 <- ((x1*wbar1)-rWD)/wbar
      x2 <- ((x2*wbar1)+rWD)/wbar
      x3 <- ((x3*wbar3)+rWD)/wbar
      x4 <- ((x4*wbar3)-rWD)/wbar
    }
    p2 <- x1 + x3
    q2 <- 1- p2
    p2_final <- c(p2_final, q2)
    x1<-ini_x1T
    x2<-ini_x2T
    x3<-ini_x3T
    x4<-ini_x4T
  }
}

# table_topSelectedlocus<-data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop)
table_topSelectedlocus<-rbind(table_topSelectedlocus, data.frame(dis, p2_final, sel=s, Ne=N, A1=A1_value, population=pop))

table_topSelectedlocus %>% 
  mutate(population=factor(population, levels=c("Top pop.", "Bottom pop."))) %>% 
  ggplot(aes(dis,p2_final, colour=A1, group=A1)) +
  geom_line() +
  ylim(c(0,1)) +
  # xlim(c(0,0.1)) +
  facet_grid(as.factor(sel) ~ population) + 
  labs(y = "Final frequency neutral locus", 
       x = "Distance to selected locus (r)", 
       colour= "Initial Fq. A1B1 Genotype") +
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_line(colour = "grey70"), 
        panel.grid.minor = element_blank(), #element_line(colour = "grey70"), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom",
        # legend.position="none", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) +
  ggsave(paste0("plots/22_7_LDsimulation_TopSel_II_V1_1741521_V2_1718756_phIII_pink.png"), width = 7, height = 6, dpi = 400)

###

##### Mean SNP frequency:

phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))



std <- function(x) sd(x)/sqrt(length(x))

phaseIII_table_total_red_mean<- phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:26) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, GV_code) %>% 
  select(GV_code, treatment, pop_ID, VarFreq) %>% 
  summarise(mean_VarFreq=mean(VarFreq), 
            sd_VarFreq=std(VarFreq)) %>% 
  ungroup() %>% 
  select(treatment, GV_code, mean_VarFreq) %>%
  spread(treatment, mean_VarFreq) 


phaseIII_table_total_red_sd<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:26) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, GV_code) %>% 
  select(GV_code, treatment, pop_ID, VarFreq) %>% 
  summarise(mean_VarFreq=mean(VarFreq), 
            sd_VarFreq=std(VarFreq)) %>% 
  ungroup() %>% 
  select(treatment, GV_code, sd_VarFreq) %>%
  spread(treatment, sd_VarFreq) 


phaseIII_table_total_red_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_red_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean)) %>% 
  ggplot(aes(Allo_B_mean, Allo_T_mean, colour=ratio_B_T)) +
  geom_point(size=2, alpha=0.8) + 
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd, 
                     colour=ratio_B_T),
                 height=0.4, alpha=0.5) +
  geom_errorbar(aes(ymin=Allo_T_mean-Allo_T_sd,
                    ymax=Allo_T_mean+Allo_T_sd, 
                    colour=ratio_B_T),
                width=0.4, alpha=0.5) + 
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(breaks = seq(0,100,20), labels = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,100,20), labels = seq(0,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray50",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  xlab("Mean SNP Fq. Allo. Bottom") + 
  ylab("Mean SNP Fq. Allo. Top") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_allo_red.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_allo_red.svg", width = 7, height = 7)



# Sympatric populations:

phaseIII_table_total_red_mean_sym<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("SYM")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:14) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(GV_code) %>% 
  select(GV_code, VarFreq) %>% 
  summarise(num_pop=n(), 
            mean_VarFreq_sym=mean(VarFreq), 
            sd_VarFreq_sym=std(VarFreq)) 

head(phaseIII_table_total_red_mean_sym)


phaseIII_table_total_red_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_red_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_red_mean_sym, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_sym=if_else(is.na(mean_VarFreq_sym), 0, mean_VarFreq_sym), 
         sd_VarFreq_sym=if_else(is.na(sd_VarFreq_sym), 0, sd_VarFreq_sym)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_sym-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_sym-Allo_T_mean) %>% 
  ggplot(aes(Allo_B_mean, mean_VarFreq_sym, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0) + 
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_sym-sd_VarFreq_sym,
                    ymax=mean_VarFreq_sym+sd_VarFreq_sym,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  xlab("Mean SNP Fq. Allo. Bottom") +
  ylab("Mean SNP Fq. Sym.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) #+
  # ggsave("plots/24_meanSVFq_alloB_sym_red.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloB_sym_red.svg", width = 7, height = 7)


phaseIII_table_total_red_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_red_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_red_mean_sym, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_sym=if_else(is.na(mean_VarFreq_sym), 0, mean_VarFreq_sym), 
         sd_VarFreq_sym=if_else(is.na(sd_VarFreq_sym), 0, sd_VarFreq_sym)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_sym-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_sym-Allo_T_mean) %>% 
  ggplot(aes(Allo_T_mean, mean_VarFreq_sym, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_T_mean-Allo_T_sd,
                     xmax=Allo_T_mean+Allo_T_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_sym-sd_VarFreq_sym,
                    ymax=mean_VarFreq_sym+sd_VarFreq_sym,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  xlab("Mean SNP Fq. Allo. Top") +
  ylab("Mean SNP Fq. Sym.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_alloT_sym_red.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloT_sym_red.svg", width = 7, height = 7)


# LM populations:

phaseIII_table_total_red_mean_lm<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("LM")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:15) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(GV_code) %>% 
  select(GV_code, VarFreq) %>% 
  summarise(num_pop=n(), 
            mean_VarFreq_lm=mean(VarFreq), 
            sd_VarFreq_lm=std(VarFreq)) 

head(phaseIII_table_total_red_mean_lm)


phaseIII_table_total_red_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_red_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_red_mean_lm, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_lm=if_else(is.na(mean_VarFreq_lm), 0, mean_VarFreq_lm), 
         sd_VarFreq_lm=if_else(is.na(sd_VarFreq_lm), 0, sd_VarFreq_lm)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_lm-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_lm-Allo_T_mean) %>% 
  ggplot(aes(Allo_B_mean, mean_VarFreq_lm, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0) + 
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_lm-sd_VarFreq_lm,
                    ymax=mean_VarFreq_lm+sd_VarFreq_lm,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  xlab("Mean SNP Fq. Allo. Bottom") +
  ylab("Mean SNP Fq. LM.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_alloB_lm_red.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloB_lm_red.svg", width = 7, height = 7)


phaseIII_table_total_red_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_red_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_red_mean_lm, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_lm=if_else(is.na(mean_VarFreq_lm), 0, mean_VarFreq_lm), 
         sd_VarFreq_lm=if_else(is.na(sd_VarFreq_lm), 0, sd_VarFreq_lm)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_lm-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_lm-Allo_T_mean) %>% 
  ggplot(aes(Allo_T_mean, mean_VarFreq_lm, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_T_mean-Allo_T_sd,
                     xmax=Allo_T_mean+Allo_T_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_lm-sd_VarFreq_lm,
                    ymax=mean_VarFreq_lm+sd_VarFreq_lm,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  xlab("Mean SNP Fq. Allo. Top") +
  ylab("Mean SNP Fq. LM.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_alloT_lm_red.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloT_lm_red.svg", width = 7, height = 7)



##### Mean SNP frequency:
# pink populations:


phaseIII_table_total_red<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3R" ) 

phaseIII_table_total_red$time_groupGV<-"phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% parental_phaseIII_red]<-"parental_phaseIII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_phaseII_red]<-"parental_phaseII"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_phaseI_red]<-"phaseI"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% polymorphic_parental_GV_red_20]<-"polymorphic_parental"
phaseIII_table_total_red$time_groupGV[phaseIII_table_total_red$GV_code %in% fix_parental_GV_red_20]<-"fix_parental"
phaseIII_table_total_red$time_groupGV <- factor(phaseIII_table_total_red$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))



std <- function(x) sd(x)/sqrt(length(x))

phaseIII_table_total_red_mean<- phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:26) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, GV_code) %>% 
  select(GV_code, treatment, pop_ID, VarFreq) %>% 
  summarise(mean_VarFreq=mean(VarFreq), 
            sd_VarFreq=std(VarFreq)) %>% 
  ungroup() %>% 
  select(treatment, GV_code, mean_VarFreq) %>%
  spread(treatment, mean_VarFreq) 


phaseIII_table_total_red_sd<-phaseIII_table_total_red %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:26) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, GV_code) %>% 
  select(GV_code, treatment, pop_ID, VarFreq) %>% 
  summarise(mean_VarFreq=mean(VarFreq), 
            sd_VarFreq=std(VarFreq)) %>% 
  ungroup() %>% 
  select(treatment, GV_code, sd_VarFreq) %>%
  spread(treatment, sd_VarFreq) 


phaseIII_table_total_red_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_red_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean)) %>% 
  ggplot(aes(Allo_B_mean, Allo_T_mean, colour=ratio_B_T)) +
  geom_point(size=2, alpha=0.8) + 
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd, 
                     colour=ratio_B_T),
                 height=0.4, alpha=0.5) +
  geom_errorbar(aes(ymin=Allo_T_mean-Allo_T_sd,
                    ymax=Allo_T_mean+Allo_T_sd, 
                    colour=ratio_B_T),
                width=0.4, alpha=0.5) +
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(breaks = seq(0,100,20), labels = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,100,20), labels = seq(0,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray50",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar") + 
  xlab("Mean SNP Fq. Allo. Bottom") + 
  ylab("Mean SNP Fq. Allo. Top") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) 
# ggsave("plots/24_meanSVFq_allo_red.png", width = 7, height = 7, dpi = 400)
# ggsave("plots/24_meanSVFq_allo_red.svg", width = 7, height = 7)





### pink populations:

phaseIII_table_total_pink<-table_total %>% 
  filter(!(sample %in% as.vector(as.numeric(c(lane_5, lane_6)))) & 
           (time!=(-16) | VarFreq>70) &
           (time!=(-10) | VarFreq>70) &
           time %in% c(-16,-10,-5,0,53) &
           parental == "P3P" ) 

phaseIII_table_total_pink$time_groupGV<-"phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% parental_phaseIII_pink]<-"parental_phaseIII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_phaseII_pink]<-"parental_phaseII"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_phaseI_pink]<-"phaseI"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% polymorphic_parental_GV_pink_20]<-"polymorphic_parental"
phaseIII_table_total_pink$time_groupGV[phaseIII_table_total_pink$GV_code %in% fix_parental_GV_pink_20]<-"fix_parental"
phaseIII_table_total_pink$time_groupGV <- factor(phaseIII_table_total_pink$time_groupGV, levels = c("fix_parental", "polymorphic_parental", "phaseI", "parental_phaseII", "parental_phaseIII", "phaseIII"))



std <- function(x) sd(x)/sqrt(length(x))

phaseIII_table_total_pink_mean<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:25) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, GV_code) %>% 
  select(GV_code, treatment, pop_ID, VarFreq) %>% 
  summarise(mean_VarFreq=mean(VarFreq), 
            sd_VarFreq=std(VarFreq)) %>% 
  ungroup() %>% 
  select(treatment, GV_code, mean_VarFreq) %>%
  spread(treatment, mean_VarFreq) 


phaseIII_table_total_pink_sd<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  # filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("Allo_T", "Allo_B")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:25) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(treatment, GV_code) %>% 
  select(GV_code, treatment, pop_ID, VarFreq) %>% 
  summarise(mean_VarFreq=mean(VarFreq), 
            sd_VarFreq=std(VarFreq)) %>% 
  ungroup() %>% 
  select(treatment, GV_code, sd_VarFreq) %>%
  spread(treatment, sd_VarFreq) 


phaseIII_table_total_pink_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_pink_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean)) %>% 
  ggplot(aes(Allo_B_mean, Allo_T_mean, colour=ratio_B_T)) +
  geom_point(size=2, alpha=0.8) + 
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd, 
                     colour=ratio_B_T),
                 height=0.4, alpha=0.5) +
  geom_errorbar(aes(ymin=Allo_T_mean-Allo_T_sd,
                    ymax=Allo_T_mean+Allo_T_sd, 
                    colour=ratio_B_T),
                width=0.4, alpha=0.5) +
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(breaks = seq(0,100,20), labels = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,100,20), labels = seq(0,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray50",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar") + 
  xlab("Mean SNP Fq. Allo. Bottom") + 
  ylab("Mean SNP Fq. Allo. Top") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_allo_pink.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_allo_pink.svg", width = 7, height = 7)


#####

# Sympatric populations:

phaseIII_table_total_pink_mean_sym<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("SYM")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:15) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(GV_code) %>% 
  select(GV_code, VarFreq) %>% 
  summarise(num_pop=n(), 
            mean_VarFreq_sym=mean(VarFreq), 
            sd_VarFreq_sym=std(VarFreq)) 

head(phaseIII_table_total_pink_mean_sym)


phaseIII_table_total_pink_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_pink_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_pink_mean_sym, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_sym=if_else(is.na(mean_VarFreq_sym), 0, mean_VarFreq_sym), 
         sd_VarFreq_sym=if_else(is.na(sd_VarFreq_sym), 0, sd_VarFreq_sym)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_sym-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_sym-Allo_T_mean) %>% 
  ggplot(aes(Allo_B_mean, mean_VarFreq_sym, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0) + 
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_sym-sd_VarFreq_sym,
                    ymax=mean_VarFreq_sym+sd_VarFreq_sym,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar") +
  xlab("Mean SNP Fq. Allo. Bottom") +
  ylab("Mean SNP Fq. Sym.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  ggsave("plots/24_meanSVFq_alloB_sym_pink.png", width = 7, height = 7, dpi = 400)
# ggsave("plots/24_meanSVFq_alloB_sym_pink.svg", width = 7, height = 7)


phaseIII_table_total_pink_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_pink_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_pink_mean_sym, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_sym=if_else(is.na(mean_VarFreq_sym), 0, mean_VarFreq_sym), 
         sd_VarFreq_sym=if_else(is.na(sd_VarFreq_sym), 0, sd_VarFreq_sym)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_sym-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_sym-Allo_T_mean) %>% 
  ggplot(aes(Allo_T_mean, mean_VarFreq_sym, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_T_mean-Allo_T_sd,
                     xmax=Allo_T_mean+Allo_T_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_sym-sd_VarFreq_sym,
                    ymax=mean_VarFreq_sym+sd_VarFreq_sym,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar") +
  xlab("Mean SNP Fq. Allo. Top") +
  ylab("Mean SNP Fq. Sym.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_alloT_sym_pink.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloT_sym_pink.svg", width = 7, height = 7)


# LM populations:

phaseIII_table_total_pink_mean_lm<-phaseIII_table_total_pink %>% 
  filter(time_groupGV!="fix_parental") %>% 
  filter(time_groupGV!="phaseIII") %>%
  # filter(time %in% c(53)) %>%
  filter(treatment %in% c("LM")) %>% 
  filter(time %in% c(53)) %>% 
  mutate(time_pop_ID=paste0(treatment, "__", time,"__",pop_treatment)) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq, Pvalue) %>% 
  group_by(Chrom, Position, GV_code, time_groupGV, time_pop_ID) %>% 
  dplyr::summarize(VarFreq=VarFreq[which(Pvalue==min(Pvalue))][1]) %>% 
  select(Chrom, Position, GV_code, time_groupGV, time_pop_ID, VarFreq) %>% 
  spread(time_pop_ID, VarFreq, fill=0) %>% 
  gather("time_pop_ID", "VarFreq", 5:15) %>% 
  separate(time_pop_ID, c("treatment", "time", "pop_ID"), "__") %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R", "SYM", "LM", "Para_T", "Para_B")) %>% 
  # filter(treatment %in% c("Allo_T", "Allo_B", "P3R")) %>% 
  group_by(GV_code) %>% 
  select(GV_code, VarFreq) %>% 
  summarise(num_pop=n(), 
            mean_VarFreq_lm=mean(VarFreq), 
            sd_VarFreq_lm=std(VarFreq)) 

head(phaseIII_table_total_pink_mean_lm)


phaseIII_table_total_pink_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_pink_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_pink_mean_lm, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_lm=if_else(is.na(mean_VarFreq_lm), 0, mean_VarFreq_lm), 
         sd_VarFreq_lm=if_else(is.na(sd_VarFreq_lm), 0, sd_VarFreq_lm)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_lm-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_lm-Allo_T_mean) %>% 
  ggplot(aes(Allo_B_mean, mean_VarFreq_lm, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0) + 
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_B_mean-Allo_B_sd,
                     xmax=Allo_B_mean+Allo_B_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_lm-sd_VarFreq_lm,
                    ymax=mean_VarFreq_lm+sd_VarFreq_lm,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar") +
  xlab("Mean SNP Fq. Allo. Bottom") +
  ylab("Mean SNP Fq. LM.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_alloB_lm_pink.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloB_lm_pink.svg", width = 7, height = 7)


phaseIII_table_total_pink_mean %>% 
  rename(Allo_B_mean=Allo_B, 
         Allo_T_mean=Allo_T) %>% 
  merge(phaseIII_table_total_pink_sd, by="GV_code") %>% 
  rename(Allo_B_sd=Allo_B, 
         Allo_T_sd=Allo_T) %>% 
  mutate(ratio_B_T=Allo_B_mean/(Allo_B_mean+Allo_T_mean), 
         diff_B_T=Allo_B_mean-Allo_T_mean) %>% 
  merge(phaseIII_table_total_pink_mean_lm, by="GV_code", all.x = T) %>% 
  mutate(mean_VarFreq_lm=if_else(is.na(mean_VarFreq_lm), 0, mean_VarFreq_lm), 
         sd_VarFreq_lm=if_else(is.na(sd_VarFreq_lm), 0, sd_VarFreq_lm)) %>% 
  mutate(diff_Sym_B=mean_VarFreq_lm-Allo_B_mean, 
         diff_Sym_T=mean_VarFreq_lm-Allo_T_mean) %>% 
  ggplot(aes(Allo_T_mean, mean_VarFreq_lm, colour=ratio_B_T)) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  geom_point(size=2, alpha=0.9) +
  geom_errorbarh(aes(xmin=Allo_T_mean-Allo_T_sd,
                     xmax=Allo_T_mean+Allo_T_sd,
                     colour=ratio_B_T),
                 height=0.4, alpha=0.7) +
  geom_errorbar(aes(ymin=mean_VarFreq_lm-sd_VarFreq_lm,
                    ymax=mean_VarFreq_lm+sd_VarFreq_lm,
                    colour=ratio_B_T),
                width=0.4, alpha=0.7) +
  # geom_abline(intercept = 0, slope = -1) + 
  scale_x_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = seq(-1,1,0.2)) +
  scale_colour_gradient2(low = "cornflowerblue", mid = "gray90",
                         high = "red3", midpoint = 0.5, space = "Lab",
                         na.value = "grey50", guide = "colourbar") +
  xlab("Mean SNP Fq. Allo. Top") +
  ylab("Mean SNP Fq. LM.") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none",
        axis.title = element_text(size=20, colour="black"), 
        axis.text.x=element_text(size=20, colour="black"), 
        axis.text.y = element_text(colour="black", size=20)) +
  # ggsave("plots/24_meanSVFq_alloT_lm_pink.png", width = 7, height = 7, dpi = 400)
  ggsave("plots/24_meanSVFq_alloT_lm_pink.svg", width = 7, height = 7)




####
