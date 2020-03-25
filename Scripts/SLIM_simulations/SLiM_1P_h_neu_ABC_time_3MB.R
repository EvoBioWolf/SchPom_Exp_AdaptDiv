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
  " ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_h_neu_ABC_time_3MB.txt"), stdout=T)
#as.numeric(output[length(output)])
#return(c(gen, as.numeric(output[length(output)-1]),as.numeric(output[length(output)])))
return(output)
}

countFq<-function(vect, lowV, highV){
	countsFq<-length(vect[vect>lowV & vect<=highV])
	return(countsFq)
}

final_table<-data.frame(sim=0, fq=0, 
  counts_neutral=0)

#for(sim in seq(1,Nsim)){
a<-runSLiM(Ngenerations)
# vector of allele frequencies for polymorphic variants:
neu_freq_vector<-as.numeric(a[(which(a=="Neutral Freq Effect size:")+1):(which(a=="EffectSubs")-1)])*2
# check number of fix variants
if (which(a=="EffectSubs") != length(a)){
  effect_subs<-as.numeric(a[(which(a=="EffectSubs")+1):length(a)])
} else {
  effect_subs<-c()
}
# Frequencies table:
if(length(effect_subs)>0){
  neu_freq_vector2<-c(neu_freq_vector,rep(1,sum(effect_subs == 0)))
} else {
  neu_freq_vector2<-c(neu_freq_vector)
}
for(fq_used in seq(0,0.995, 0.005)){
  final_table<-rbind(final_table, c(sim=Nsim, fq=fq_used, 
    counts_neutral=countFq(neu_freq_vector2,fq_used, fq_used+0.005)))
}
final_table<-final_table[-1,]
write.table(final_table, file = paste0(outputFile,"_fq.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
# summary table:
final_table_summary<-data.frame(sim=Nsim, 
P_mutations_neutral=length(neu_freq_vector), 
substitutions=length(effect_subs))
write.table(final_table_summary, file = paste0(outputFile,"_summary.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
# empty tables for next simulation:
final_table<-data.frame(sim=0, fq=0, 
counts_neutral=0, 
counts_sel=0)
#}

