library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
setwd('~/Desktop/_mtree/_simulated')
source('seqsim.R')
load('RV217bestMCCtrees_genes.RData')

#MODELTEST
k80<-lapply(1:10,function(l){
  print(l)
  seqsim(tre,length=20,model='K80',alpha=10,beta=1,basefreq=c(0.4,0.3,0.2,0.1),option='normal')
})

hky<-lapply(1:10,function(l){
  print(l)
  seqsim(tre,length=20,model='HKY',alpha=10,beta=1,basefreq=c(0.4,0.3,0.2,0.1),option='normal')
})

gtr<-lapply(1:10,function(l){
  print(l)
  seqsim(tre,length=20,model='GTR',basefreq=c(0.4,0.3,0.2,0.1),option='normal')
})


mods<-c(k80,hky,gtr)
mod<-c(rep("k80",10),rep("hky",10),rep("gtr",10))


lapply(1:length(mods),function(i){
  saveAlignment(mods[[i]],file=paste(mod[i],i,".fas",sep=""),skip.internal=T)
})		