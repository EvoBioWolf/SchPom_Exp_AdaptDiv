#!/usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)

# parameters for a test run:
Nsim<-1
Ngenerations<-750
K<-3e5
fracSel<-3e3
L<-1e6
H<-0.01
selection<-0.1
propSused<-1000
mut_rate<-2e-8
recon_rate<-1e-4
outputFile<-"test"

# parameters taken from the batch run:
Nsim<-args[1]
Ngenerations<-args[2]
selection<-args[3]
propSused<-args[4]
L<-as.numeric(args[5])*1e6
outputFile<-args[6]


runSLiM <- function(x) {
gen <- x
#cat("Running SLiM with seed ", seed, ", mu = ", mu, "\n");
output <- system2("/proj/uppstore2017159/b2014286/private/tools/SLiM/build/slim", c("-d", paste0("gen=", gen), 
  "-d", paste0("K=", K), 
  "-d", paste0("fracSel=", fracSel), 
  "-d", paste0("L=", L), 
  "-d", paste0("H=", H), 
  "-d", paste0("sel=", selection), 
  "-d", paste0("propS=", propSused), 
  "-d", paste0("mut_rate=", mut_rate), 
  "-d", paste0("recon_rate=", recon_rate), 
  " ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_h_sel_ABC_time_3MB.txt"), stdout=T)
#as.numeric(output[length(output)])
#return(c(gen, as.numeric(output[length(output)-1]),as.numeric(output[length(output)])))
return(output)
}

countFq<-function(vect, lowV, highV){
	countsFq<-length(vect[vect>lowV & vect<=highV & !(is.na(vect))])
	return(countsFq)
}

final_table<-data.frame(sim=0, fq=0, 
  counts_neutral=0, 
  counts_sel=0)

final_table_effect<-data.frame(sim=0, fq=0, 
  counts_effect=0)

#for(sim in seq(1,Nsim)){
a<-runSLiM(Ngenerations)
# P_mutations=as.numeric(a[16])
# P_mutations_neutral=as.numeric(a[18])
# P_mutations_sel=as.numeric(a[20])
# substitutions=as.numeric(a[22])
# substitutions_neutral=as.numeric(a[24])
# substitutions_sel=as.numeric(a[26])
neu_freq_vector<-as.numeric(a[(which(a=="Neutral Freq Effect size:")+1):(which(a=="EffectSubs")-1)])*2
#neu_effect_vector<-as.numeric(a[(which(a=="Neutral Effect size:")+1):(which(a=="Neutral Freq Effect size:")-1)])
#if (length(a[(which(a=="Sel Freq Effect size:")+1):(which(a=="Neutral Effect size:")-1)])>0){
if ((which(a=="Sel Freq Effect size:")+1)!=(which(a=="Neutral Effect size:"))) {
  sel_freq_vector<-as.numeric(a[(which(a=="Sel Freq Effect size:")+1):(which(a=="Neutral Effect size:")-1)])*2
} else {
  sel_freq_vector<-c()
}
if((which(a=="Sel Effect size:")+1)!=(which(a=="Sel Freq Effect size:"))){
  sel_effect_vector<-as.numeric(a[(which(a=="Sel Effect size:")+1):(which(a=="Sel Freq Effect size:")-1)])
} else {
  sel_effect_vector<-c()
}
if (which(a=="EffectSubs") != length(a)){
  effect_subs<-as.numeric(a[(which(a=="EffectSubs")+1):length(a)])
} else {
  effect_subs<-c()
}
# Frequencies table:
if(length(effect_subs)>0 & sum(effect_subs == 0)>0){
    neu_freq_vector2<-c(neu_freq_vector,rep(1,sum(effect_subs == 0)))
} else {
  neu_freq_vector2<-c(neu_freq_vector)
}
if(length(effect_subs)>0 & sum(effect_subs != 0)>0){
  sel_freq_vector2<-c(sel_freq_vector,rep(1,sum(effect_subs != 0)))
} else {
  sel_freq_vector2<-c(sel_freq_vector)
}
for(fq_used in seq(0,0.995, 0.005)){
  final_table<-rbind(final_table, c(sim=Nsim, fq=fq_used, 
    counts_neutral=countFq(neu_freq_vector2,fq_used, fq_used+0.005), 
    counts_sel=countFq(sel_freq_vector2,fq_used, fq_used+0.005)))
}
final_table<-final_table[-1,]
write.table(final_table, file = paste0(outputFile,"_fq.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
# Effect table
if(length(effect_subs)>0){
  sel_effect_vector2<-c(sel_effect_vector,effect_subs)
} else {
  sel_effect_vector2<-c(sel_effect_vector)
}
for(fq_used in seq(0,0.995, 0.005)){
  final_table_effect<-rbind(final_table_effect, c(sim=Nsim, fq=fq_used, 
    counts_effect=countFq(sel_effect_vector2,fq_used, fq_used+0.005))) 
}
final_table_effect<-final_table_effect[-1,]
write.table(final_table_effect, file = paste0(outputFile,"_effect.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
# summary table:
final_table_summary<-data.frame(sim=Nsim, 
P_mutations=length(c(sel_freq_vector, neu_freq_vector)), 
P_mutations_neutral=length(neu_freq_vector), 
P_mutations_sel=length(sel_freq_vector), 
substitutions=length(effect_subs), 
substitutions_neutral=sum(effect_subs == 0), 
substitutions_sel=sum(effect_subs != 0))
write.table(final_table_summary, file = paste0(outputFile,"_summary.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
# empty tables for next simulation:
final_table<-data.frame(sim=0, fq=0, 
counts_neutral=0, 
counts_sel=0)
final_table_effect<-data.frame(sim=0, fq=0, 
counts_effect=0)
#}

