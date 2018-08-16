library(RPANDA)
library(igraph)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
setwd('~/Desktop/_mtree/_simulated/_nucsim')
source('seqsim.R')
source('~/Documents/GitHub/pflam/_R/_JSD.R')
setwd('~/Desktop/_mtree/_simulated/_nucsim/_titvtest')
tr<-read.tree('gtrTree')
#load('RV217bestMCCtrees_genes.RData')
#load('genetree_for_simulation.RData')

#RATETEST	
##nucleotide
alpha<-c(50,25,1,1,50)#transition
beta<-c(1,25,50,1,50)#transversion
	c()->hky
	for(i in 1:length(alpha)){
		print(i)
hky[[i]]<-lapply(1:10,function(l){
		print(l)
	seqsim(tr,length=666,model='HKY',option="normal",
		basefreq=rep(0.25,4),alpha=alpha[i],beta=beta[i])->t
	saveAlignment(t,file=paste("HKY",i,l,".fas",sep=""),skip.internal=T)
	rm(t)	
	})
}

####################ANALYSIS_OF_RATETREES#####################################ANALYSIS_OF_RATETREES#####################################ANALYSIS_OF_RATETREES#####################################ANALYSIS_OF_RATETREES########
library(igraph)
source('~/Desktop/_mtree/_R/_JSD.R')
files<-list.files(path="~/Desktop/_mtree/_simulated/_ratetest/_ratetestrees",pattern="*.treefile")
	setwd('~/Desktop/_mtree/_simulated/_ratetest/_ratetestrees')
	ratetree<-lapply(files,read.tree)

alignments<-list.files(path="~/Desktop/_mtree/_simulated/_ratetest/_ratetestrees",pattern="*.fasta")
	ratefas<-lapply(alignments,read.alignment,format='fasta')
	sel<-lapply(c(7,9,10),function(q){
		print(q)
		kaks(ratefas[[q]])->r
		c(mean(r$ka/r$ks),min(r$ka/r$ks),max(r$ka/r$ks))
	})

##compute spectral density profile
#set SDP function
gete<-function(phy){
eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=T),
	only.values=T)$values->x
	subset(x,x>1)->e
	return(e)
}

e<-lapply(ratetree,gete)

c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(e)){
	dens(log(e[[n]]))->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}

alpha<-c(rep('alpha10',10),rep('alpha5',10),rep('alpha1',10))
beta<-c(rep('beta1',10),rep('beta5',10),rep('beta10',10))
tab<-cbind(log(pe),skew,height,alpha,beta)
	colnames(tab)<-c("pe","skew","height","alpha","beta")
write.table(tab,file="ratetest_table")


Ds<-c() 
for(i in 1:length(e)){
     Ds<-as.data.frame(cbind(Ds,abs(d[[i]]$x)))
     }
virus_JSD <- as.matrix(JSDist(Ds))
	rownames(virus_JSD)<-alpha
	colnames(virus_JSD)<-beta
heatmap.2(virus_JSD,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.6,cexCol=0.6)


rat<-read.table("ratetest_table",header=T)
high<-subset(rat,rat$alpha=='alpha10')
mid<-subset(rat,rat$alpha=='alpha5')
low<-subset(rat,rat$alpha=='alpha1')


sactter.grid(rat$pe,rat$height,rat$skew,angle=130,xlab="PE",ylab="height",zlab="skew",pch=20,color=as.numeric(rat$alpha))
par(mfrow=c(1,3))
boxplot(high$pe,mid$pe,low$pe,outline=F)
boxplot(high$skew,mid$skew,low$skew,outline=F)
boxplot(high$height,mid$height,low$height,outline=F)



