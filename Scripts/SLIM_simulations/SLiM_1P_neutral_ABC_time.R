#!/usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)

Nsim<-args[1]
Ngenerations<-args[2]
outputFile<-args[3]

runSLiM <- function(x) {
gen <- x
#cat("Running SLiM with seed ", seed, ", mu = ", mu, "\n");
output <- system2("/proj/uppstore2017159/b2014286/private/tools/SLiM/build/slim", c("-d", paste0("gen=", gen), " ~/private/Uppsala/Analyses/04_Genomic_analyses/11_phaseIII/scripts/SLiM_1P_neutral_ABC_time.txt"), stdout=T)
#as.numeric(output[length(output)])
#return(c(gen, as.numeric(output[length(output)-1]),as.numeric(output[length(output)])))
return(output)
}

countFq<-function(vect, lowV, highV){
	countsFq<-length(vect[vect>lowV & vect<=highV])
	return(countsFq)
}

final_table<-data.frame(sim=0, fq=0, counts=0)

for(sim in seq(1,Nsim)){
  a<-runSLiM(Ngenerations)
  FV<-as.numeric(a[14])
  SS<-as.numeric(a[15:(length(a)-1)])
  if(FV>0){
  	SS<-c(SS,rep(1,FV))
  }
  for(fq_used in seq(0,0.995, 0.005)){
  	final_table<-rbind(final_table, c(sim=sim, fq=fq_used, counts=countFq(SS,fq_used, fq_used+0.005)))
  }
  final_table<-final_table[-1,]
  write.table(final_table, file = outputFile, append = T, quote = F, sep = "\t", row.names = F, col.names = F)
  final_table<-data.frame(sim=0, fq=0, counts=0)
}

