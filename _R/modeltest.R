library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
setwd('~/Desktop/_mtree/_simulated')
source('seqsim.R')
rcoal(30)->tre

#MODELTEST
k80<-lapply(1:10,function(l){
  print(l)
  seqsim(tre,length=20,model='K80',alpha=10,beta=1,basefreq=c(0.4,0.3,0.2,0.1),option='normal')
  saveAlignment(t,file=paste("k80",l,".fas",sep=""),skip.internal=T)
})

hky<-lapply(1:10,function(l){
  print(l)
  seqsim(tre,length=20,model='HKY',alpha=10,beta=1,basefreq=c(0.4,0.3,0.2,0.1),option='normal')
  saveAlignment(t,file=paste("hky",l,".fas",sep=""),skip.internal=T)
	rm(t)
})

gtr<-lapply(1:10,function(l){
  print(l)
  seqsim(tre,length=20,model='GTR',basefreq=c(0.4,0.3,0.2,0.1),option='normal')
  saveAlignment(t,file=paste("gtr",l,".fas",sep=""),skip.internal=T)
})

