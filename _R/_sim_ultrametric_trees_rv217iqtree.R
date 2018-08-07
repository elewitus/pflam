library(RPANDA)
library(phytools)
library(igraph)
library(gplots)
library(RColorBrewer)
library(fpc)
source('~/Documents/GitHub/pflam/_R/_JSD.R')

setwd('~/Desktop/_simulated/_proc/_nonultrametric')
cols<-c(1,4,3,2,5)*15

rtab<-read.csv('RV217neutralizer_table.csv',header=T,row.names=1)
tab<-read.table('sim_viral_phylos_table_nonultrametric.txt')

pdf('simulated_nonultrametric_trees_phylospace.pdf')
sactter.grid(log(tab$pe),tab$skew,tab$height,color=colors(1)[as.numeric(tab$proc)*15],pch=20,angle=130)
	legend("topleft",levels(tab$proc),col=colors(1)[cols],pch=20,bty="n")
dev.off()

pdf('simulated_nonultrametric_trees_pairwise.pdf',6,12)
par(mfrow=c(3,1))
plot(tab$pe,tab$height,col=colors(1)[as.numeric(tab$proc)*15],pch=20)
plot(tab$pe,tab$skew,col=colors(1)[as.numeric(tab$proc)*15],pch=20)
plot(tab$skew,tab$height,col=colors(1)[as.numeric(tab$proc)*15],pch=20)
dev.off()

subs<-lapply(1:length(levels(tab$proc)),function(j){
	subset(tab,tab$proc==levels(tab$proc)[j])
})

pdf('simulated_ultrametric_trees_phylostats.pdf',6,12)
par(mfrow=c(3,1))
c()->ps
pes<-lapply(1:length(subs),function(p){c(ps,subs[[p]]$pe)->ps})
	boxplot(pes,names=levels(tab$proc),las=2,outline=F,col=colors(1)[cols])
c()->sk
skews<-lapply(1:length(subs),function(p){c(sk,subs[[p]]$skew)->sk})
	boxplot(skews,names=levels(tab$proc),las=2,outline=F,col=colors(1)[cols])	
c()->hs
heights <-lapply(1:length(subs),function(p){c(hs,subs[[p]]$height)->hs})
	boxplot(heights,names=levels(tab$proc),las=2,outline=F,col=colors(1)[cols])
dev.off()	