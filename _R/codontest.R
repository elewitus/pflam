R

library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
library(igraph)
library(seqinr)
source('~/Desktop/_mtree/_simulated/seqsim.R')
source('~/Desktop/_mtree/_R/_JSD.R')

setwd('~/Desktop/_mtree/_simulated/_ratetest/_codon')

tr<-rtree(30)
k<-rep(2^seq(-5,4,3),10)
test<-lapply(1:length(k),function(o){
	print(o)
	codsim(tr,500,k[o])
})

lapply(1:length(test),function(a){
	saveAlignment(test[[a]],file=paste("codontest",a,".fas",sep=""),skip.internal=T)
	})

cd ./_codon

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

R
library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
library(igraph)
library(seqinr)
source('~/Desktop/_mtree/_simulated/seqsim.R')
source('~/Desktop/_mtree/_R/_JSD.R')

setwd('~/Desktop/_mtree/_simulated/_ratetest/_codon')

files<-list.files(path="~/Desktop/_mtree/_simulated/_ratetest/_codon",pattern='fasta')
al<-lapply(files[seq(1,79,2)],read.alignment,format='fasta')
dnds<-lapply(al,kaks)

rate<-sapply(1:length(dnds),function(j){
	mean(dnds[[j]]$ka/dnds[[j]]$ks,na.rm=T)
})

trf<-list.files(path="~/Desktop/_mtree/_simulated/_ratetest/_codon",pattern='treefile')
trees<-lapply(trf,read.tree)

gete<-function(phy){
eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=F),
	only.values=T)$values->x
	subset(x,x>1)->e
	return(e)
}

e<-lapply(trees,gete)

c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(e)){
	dens(log(e[[n]]))->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}

kappa<-rep(seq(0.1,4,1),10)
tab<-cbind(pe,skew,height,kappa)
write.table(tab,file="Simulated_codonKappa_MGLtable.txt")
########################################################################################################################################################################################################################################################################################################################################
setwd('~/Desktop/_mtree/_simulated/_ratetest/_codon')
tab<-read.table('Simulated_codonKappa_MGLtable.txt',header=T)

sactter.grid(tab$pe,tab$skew,tab$height,pch=20,color=round(tab$kappa)+1)

kap<-lapply(1:4,function(j){
	subset(tab,tab$kappa==seq(0.1,4,1)[j])
})



	