library(RPANDA)
library(phytools)
library(igraph)
library(gplots)
library(RColorBrewer)
library(fpc)
source('~/Documents/GitHub/pflam/_R/_JSD.R')

setwd('~/Desktop/_simulated/_proc/_ultrametric')
cols<-c(1,5,4,3,2,6)*15

btab<-read.table('RV217_AlignedByParticipant_genesBeast_MGLtable.txt',header=T,row.names=1)
tab<-read.table('sim_viral_phylos_table_ultrametric_large.txt')
#env only
#etab<-subset(btab,btab$proc=='env')
#tabs<-rbind(tab,etab[,c(1,3,2,5)])
#all genes
etab<-cbind(btab[,c(1,3,2)],rep('rv217',dim(btab)[1]))
colnames(etab)<-colnames(tab)
tabs<-rbind(tab,etab)
factor(tabs$proc)->tabs$proc

pdf('simulated_ultrametric_trees_RV217beast_phylospace.pdf')
sactter.grid(tabs$pe,log(tabs$skew),log(tabs$height),color=colors(1)[as.numeric(tabs$proc)*15],pch=c(rep(1,1000),rep(20,387)),angle=230)
	legend("topleft",levels(tabs$proc),col=colors(1)[cols],pch=20,bty="n")
dev.off()

pdf('simulated_ultrametric_trees_RV217beast_pairwise.pdf',6,12)
par(mfrow=c(3,1))
plot(tabs$pe,log(tabs$height),col=colors(1)[as.numeric(tabs$proc)*15],pch=20)
plot(tabs$pe,log(tabs$skew),col=colors(1)[as.numeric(tabs$proc)*15],pch=20)
plot(log(tabs$skew),log(tabs$height),col=colors(1)[as.numeric(tabs$proc)*15],pch=20)
dev.off()

subs<-lapply(1:length(levels(tabs$proc)),function(j){
	subset(tabs,tabs$proc==levels(tabs$proc)[j])
})

pdf('simulated_ultrametric_trees_RV217beast_phylostats.pdf',6,12)
par(mfrow=c(3,1))
c()->ps
pes<-lapply(1:length(subs),function(p){c(ps,subs[[p]]$pe)->ps})
	boxplot(pes,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])
c()->sk
skews<-lapply(1:length(subs),function(p){c(sk,subs[[p]]$skew)->sk})
	boxplot(skews,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])	
c()->hs
heights <-lapply(1:length(subs),function(p){c(hs,subs[[p]]$height)->hs})
	boxplot(heights,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])
dev.off()	