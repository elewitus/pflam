library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
setwd('~/Desktop/_mtree/_simulated')
source('seqsim.R')
load('RV217bestMCCtrees_genes.RData')

#RATETEST	
setwd('~/Desktop/_mtree/_simulated/_ratetest')
alpha<-c(10,5,1)
beta<-c(1,5,10)
	c()->hky
	for(i in 1:length(alpha)){
hky[[i]]<-lapply(1:10,function(l){
	seqsim(genetrees[[l]],length=30,model='HKY',option="normal",
		basefreq=rep(0.25,4),alpha=alpha[i],beta=beta[i])
	})
}

hkys<-c(hky[[1]],hky[[2]],hky[[3]])
	
lapply(1:length(hkys),function(i){
    saveAlignment(hkys[[i]],file=paste("HKYs",i,".fas",sep=""),skip.internal=T)
  })
