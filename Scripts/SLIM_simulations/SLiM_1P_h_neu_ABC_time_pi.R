#!/usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)

# parameters for a test run:
# Nsim<-1
# Ngenerations<-750
# K<-3e6
# fracSel<-3e4
# L<-1e6
# H<-0.01
# mut_rate<-2e-9
# recon_rate<-1e-5
# outputFile<-"test"

# parameters for a test run:
Nsim<-1
Ngenerations<-750
K<-3e5
fracSel<-3e3
L<-1e6
H<-0.01
mut_rate<-2e-8
recon_rate<-1e-4
outputFile<-"test"


# parameters taken from the batch run:
Nsim<-args[1]
Ngenerations<-args[2]
L<-as.numeric(args[3])*1e6
outputFile<-args[4]
K<-args[5]
fracSel<-args[6]


runSLiM <- function(x) {
gen <- x
#cat("Running SLiM with seed ", seed, ", mu = ", mu, "\n");
output <- system2("/proj/uppstore2017159/b2014286/private/tools/SLiM/build/slim", c("-d", paste0("gen=", gen), 
  "-d", paste0("K=", K), 
  "-d", paste0("fracSel=", fracSel), 
  "-d", paste0("L=", L), 
  "-d", paste0("H=", H), 
  "-d", paste0("mut_rate=", mut_rate), 
  "-d", paste0("recon_rate=", recon_rate), 
  " ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_h_neu_ABC_time_pi.txt"), stdout=T)
#as.numeric(output[length(output)])
#return(c(gen, as.numeric(output[length(output)-1]),as.numeric(output[length(output)])))
return(output)
}


#for(sim in seq(1,Nsim)){
a<-runSLiM(Ngenerations)
# vector of allele frequencies for polymorphic variants:
neu_freq_vector<-as.numeric(a[(which(a=="Neutral Freq Effect size:")+1):(which(a=="EffectSubs")-1)])*2
neu_freq_vector2<-neu_freq_vector[neu_freq_vector>0.05]
pi_value<-mean(((neu_freq_vector*100)*(100-(neu_freq_vector*100)))/(((100-1)*100)/2))
pi_value2<-mean(((neu_freq_vector2*100)*(100-(neu_freq_vector2*100)))/(((100-1)*100)/2))

# check number of fix variants
final_table_summary<-data.frame(sim=Nsim, 
  K=K,
  pi_mean=pi_value, 
  pi5_mean=ifelse(is.na(pi_value2), 0, pi_value2))

write.table(final_table_summary, file = paste0(outputFile,"__pi.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)


