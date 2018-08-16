#R

library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
library(igraph)
library(seqinr)
library(gtools)
source('~/Desktop/_mtree/_simulated/_nucsim/seqsim.R')
source('~/Desktop/_mtree/_R/_JSD.R')

setwd('~/Desktop/_mtree/_simulated/_nucsim/_ratetest/_dnds')



tr<-rcoal(20)

##dNdS variation
#dnds<-rep(2^seq(-6,3,3),2)
dnds<-rep(seq(0.1,2,0.8),10)
lapply(27:length(dnds),function(o){
	print(o)
	codsim(tr,300,1,0.9*dnds[o],dnds[o],1,plot=F)->t
	saveAlignment(t,file=paste("dndstest",o,".fas",sep=""),skip.internal=T)
	rm(t)
})



##TiTv variation
#titv<-rep(2^seq(-5,4,3),10)
#test<-lapply(1:length(k),function(o){
#	print(o)
#	codsim(tr,500,titv[o],1,plot=F)
#})

#lapply(1:length(test),function(a){
#	saveAlignment(test[[a]],file=paste("titvtest",a,".fas",sep=""),skip.internal=T)
#	})



##align with mafft
for fasta_file in $(ls *.fas)
do
linsi $fasta_file > $fasta_file.fasta
done


##build trees
for msa_file in $(ls *.fasta)
do
./iqtree -s $msa_file
done

#files<-list.files(path="~/Desktop/_mtree/_simulated/_nucsim/_ratetest/_dnds",pattern='.treefile')
 # files<-files[seq(1,112,9)]
#al<-lapply(files,read.alignment,format='fasta')
#dnds<-lapply(al,kaks)

#rate<-sapply(c(1:5,7:11,13,14),function(j){
 # print(j)
  #kaks(al[[j]])->dnds[[j]]
	#mean(dnds[[j]]$ka/dnds[[j]]$ks,na.rm=T)
#})

trf<-list.files(path="~/Desktop/_mtree/_simulated/_nucsim/_ratetest/_dnds",pattern='.treefile')
  #trf<-trf[c(1:5,7:11,13,14)]
  mixedsort(trf)->trf
trees<-lapply(trf,read.tree)
#trees<-lapply(trees,drop.tip,c(1,2))

rescales<-function(tree,scale){
  tree$edge.length<-
    tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
}

gete<-function(phy){
abs(eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=F),
	only.values=T)$values)
}
#trees<-lapply(trees,rescales,10)
#tr<-lapply(1:length(trees),function(c){chronos(trees[[c]],0.1,model="relaxed")})
e<-lapply(trees,gete)

c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(trf)){
	dens(e[[n]])->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}


tab<-cbind(pe,skew,height,dnds)
colnames(tab)<-c('pe','skew','height','dnds')
#tab[order(tab$dnds),]->tabs
write.table(tab,file="Simulated_dnds_nostem_MGLtable.txt")
########################################################################################################################################################################################################################################################################################################################################
setwd('~/Desktop/_mtree/_simulated/_nucsim/_ratetest/_dnds')
tab<-read.table('Simulated_dnds_nostem_MGLtable.txt',header=T,row.names=1)
tab<-read.table('Simulated_dnds_long_MGLtable.txt',header=T,row.names=1)
cols<-colors(1)[c(rep(10,100),rep(20,100),rep(30,100))]


sactter.grid(tab$pe,tab$skew,tab$height,angle=230,pch=20,color=cols)
plot(tab$dnds,tab$pe,col=cols,pch=20)

rate<-c(0.1,0.7,1.9)
kap<-lapply(1:3,function(j){
	subset(tab,tab$dnds==rate[j])
})

pes<-list(kap[[1]]$pe,kap[[2]]$pe,kap[[3]]$pe)
skews<-list(kap[[1]]$skew,kap[[2]]$skew,kap[[3]]$skew)
heights<-list(kap[[1]]$height,kap[[2]]$height,kap[[3]]$height)

pdf('Simulated_dnds_nostem_phylostats_boxplot.pdf',6,12)
par(mfrow=c(3,1),mar=c(4,4,1,1))
boxplot(pes,ylab=expression(lambda),axes=F);axis(2,las=2)
stripchart(pes,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=colors(1)[seq(10,50,10)])
boxplot(skews,ylab=expression(psi),axes=F);axis(2,las=2)
stripchart(skews,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=colors(1)[seq(10,50,10)])
boxplot(heights,ylab=expression(eta),xlab='dN/dS',axes=F);axis(2,las=2)
	axis(1,at=1:4,label=seq(0.1,1.9,0.6))
stripchart(heights,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=colors(1)[seq(10,50,10)]) 
dev.off()

tr<-lapply(1:3,function(c){chronos(trees[[c]],0.1,model="relaxed")})

pdf('Simulated_dnds_nostem_phylostats_boxplot.pdf',6,12)
par(mfrow=c(3,2),mar=c(4,3,1,1))
lapply(1:3,function(p){
	drop.tip(trees[[p]],1)->tree
	multi2di(tree)->trs
	plot(ladderize(trees[[p]]),show.tip.label=F,edge.color=cols[seq(1,300,100)][p])
	plot(d[[p]]$x,dsc[[p]],type='l',ylim=c(0,1.2),xlim=c(-5,20),
		col=cols[seq(1,300,100)][p])
		polygon(d[[p]]$x,dsc[[p]],col=cols[seq(1,300,100)][p])
	})
dev.off()

bfiles<-list.files(path="~/Desktop/_mtree/_simulated/_nucsim/_ratetest/_dnds",pattern='.mcc.trees')
btrees<-lapply(bfiles,read.nexus)
#btrees<-lapply(btrees,bind.tree,tip)
btrees<-lapply(btrees,rescales,1)

bspec<-lapply(btrees,gete)

c()->bd;c()->bdsc;c()->bpe;c()->bskew;c()->bheight
for(n in 1:length(btrees)){
	dens(bspec[[n]])->bd[[n]]
	bd[[n]]$y/integr(bd[[n]]$x,bd[[n]]$y)->bdsc[[n]]
	max(bspec[[n]])->bpe[[n]]
	skewness(bdsc[[n]])->bskew[[n]]
	max(bdsc[[n]])->bheight[[n]]
	}

pdf('Simulated_dnds_trees_profiles.pdf',6,12)
par(mfrow=c(4,2))
lapply(1:length(btrees),function(b){
	cols=colors(1)[seq(10,50,10)][b]
	plot(ladderize(btrees[[b]]),show.tip.label=F,edge.color=cols,root.edge=T)
	plot(bd[[b]]$x,bdsc[[b]],type='l',col=cols,ylim=c(0,0.05))
})
dev.off()

############dNdS_inverse_relation_to_eta##########################################################################################################dNdS_inverse_relation_to_eta##################################################dNdS_inverse_relation_to_eta#####################################################################################################################################################
#define constant, linear, and exponential speciation functions
l1<-function(t,y){y[1]}
l3<-function(t,y){y[1]*exp(y[2]*t)}

#define constant, linear, and exponential extinction functions
m0<-function(t,y){0}

#define speciation parameters
lp1<-c(0.01)
lp6<-c(0.01,0.2)
lp7<-c(0.01,-0.1)

#define extinction parameters
mp0<-c(0)

e<-rep(seq(0.05,0.11,0.02),100)
utrees<-lapply(1:400,function(q){
	print(q)
	tess.sim.taxa.age(1,20,30,l3(15,c(0.07,e[q])),0,samplingProbability=0.2)[[1]]
	})				

getit<-function(phy){
eigen(
graph.laplacian(
graph.adjacency(
data.matrix(dist.nodes(phy)),
weighted=T),
normalized=F),
only.values=T)$values->x
subset(x,x>1)
}

uspec<-lapply(utrees,getit)

c()->ud;c()->udsc;c()->upe;c()->uskew;c()->uheight
for(n in 1:length(utrees)){
	dens(log(uspec[[n]]))->ud[[n]]
	ud[[n]]$y/integr(ud[[n]]$x,ud[[n]]$y)->udsc[[n]]
	max(uspec[[n]])->upe[[n]]
	skewness(udsc[[n]])->uskew[[n]]
	max(udsc[[n]])->uheight[[n]]
	}

cbind(upe,uskew,uheight,e)->tabs
colnames(tabs)<-c('pe','skew','height','rate')
write.table(tabs,file='Simulated_dnds_increaseRates_MGLtable.txt')

pdf('Simulated_dnds_trees_profiles.pdf',6,12)
par(mfrow=c(4,2))
lapply(1:4,function(b){
	cols=colors(1)[seq(10,50,10)][b]
	plot(ladderize(utrees[[b]]),show.tip.label=F,edge.color=cols,root.edge=T)
	plot(ud[[b]]$x,udsc[[b]],type='l',col=cols)#,xlim=c(-20,2500),ylim=c(0,0.002))
})
dev.off()

tab<-read.table('Simulated_dnds_increaseRates_MGLtable.txt',header=T,row.names=1)

kap<-lapply(1:4,function(j){
	subset(tab,tab$rate==tab$rate[j])
})

pes<-list(kap[[1]]$pe,kap[[2]]$pe,kap[[3]]$pe,kap[[4]]$pe)
skews<-list(kap[[1]]$skew,kap[[2]]$skew,kap[[3]]$skew,kap[[4]]$skew)
heights<-list(kap[[1]]$height,kap[[2]]$height,kap[[3]]$height,kap[[4]]$height)


pdf('Simulated_dnds_increaseRates_phylostats_boxplot.pdf',6,12)
par(mfrow=c(3,1),mar=c(4,4,1,1))
boxplot(pes,ylab=expression(lambda),axes=F);axis(2,las=2)
stripchart(pes,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=colors(1)[seq(10,50,10)])
boxplot(skews,ylab=expression(psi),axes=F);axis(2,las=2)
stripchart(skews,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=colors(1)[seq(10,50,10)])
boxplot(heights,ylab=expression(eta),xlab='dN/dS',axes=F);axis(2,las=2)
	axis(1,at=1:4,label=seq(0.1,1.9,0.6))
stripchart(heights,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=colors(1)[seq(10,50,10)]) 
dev.off()

