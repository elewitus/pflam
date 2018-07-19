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

##dNdS variation
#dnds<-rep(2^seq(-6,3,3),2)
dnds<-rep(seq(0.1,2,0.3),2)
lapply(1:length(dnds),function(o){
	print(o)
	codsim(tr,100,1,0.9*dnds[o],dnds[o],1,plot=F)->t
	saveAlignment(t,file=paste("dndstest",o,".fas",sep=""),skip.internal=T)
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
  files<-files[seq(1,112,8)]
al<-lapply(files,read.alignment,format='fasta')
#dnds<-lapply(al,kaks)

rate<-sapply(c(1:5,7:11,13,14),function(j){
  print(j)
  kaks(al[[j]])->dnds[[j]]
	mean(dnds[[j]]$ka/dnds[[j]]$ks,na.rm=T)
})

trf<-list.files(path="~/Desktop/_mtree/_simulated/_ratetest/_codon",pattern='treefile')
  trf<-trf[c(1:5,7:11,13,14)]
trees<-lapply(trf,read.tree)


rescales<-function(tree,scale){
  tree$edge.length<-
    tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
}

retrees<-lapply(trees,rescales,1)

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

e<-lapply(mcct,gete)

c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(e)){
	dens(log(e[[n]]))->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}


tab<-cbind(pe,skew,height,dnds)
write.table(tab,file="Simulated_codonKappa_MGLtable.txt")
########################################################################################################################################################################################################################################################################################################################################
setwd('~/Desktop/_mtree/_simulated/_ratetest/_codon')
tab<-read.table('Simulated_codonKappa_MGLtable.txt',header=T)

sactter.grid(tab$pe,tab$skew,tab$height,pch=20,color=round(tab$kappa)+1)

kap<-lapply(1:4,function(j){
	subset(tab,tab$kappa==seq(0.1,4,1)[j])
})

lapply(c(6:9,11,12),function(s){plot(density(e[[s]]))})
lapply(c(1:5,10),function(s){plot(density(e[[s]]))})

setwd('~/Desktop/_t')
files<-list.files(path='~/Desktop/_t',pattern='treefile')
align<-lapply(files,read.alignment,format='fasta')
tree<-lapply(files,read.tree)

c()->dd;c()->dnds
for(j in 1:length(align)){
	kaks(align[[j]])->dd[[j]]
	median(dd[[j]]$ka/dd[[j]]$ks,na.rm=T)->dnds[[j]]
	}

############dNdS_inverse_relation_to_eta##########################################################################################################dNdS_inverse_relation_to_eta##################################################dNdS_inverse_relation_to_eta#####################################################################################################################################################

sp<-lapply(tree,spectR)	
sps<-sapply(1:length(sp),function(q){
	c(max(sp[[q]]$eigenvalues),skewness(dens(sp[[q]]$eigenvalues)$x),max(dens(sp[[q]]$eigenvalues)$y))	
})

	
