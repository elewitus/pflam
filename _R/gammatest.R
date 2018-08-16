library(RPANDA)
library(igraph)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
setwd('~/Desktop/_mtree/_simulated/_nucsim')
source('seqsim.R')
source('~/Documents/GitHub/pflam/_R/_JSD.R')
load('RV217bestMCCtrees_genes.RData')
setwd('~/Desktop/_mtree/_simulated/_nucsim/_gammatest')
tr<-read.tree('gtrTree')
#genetree<-rep(c(genetrees[[1]],genetrees[[3]]),50)

#################simulate_alignments_under_different_gamma_profiles##################################simulate_alignments_under_different_gamma_profiles##

gtrNormal<-lapply(1:20,function(l){
	print(l)
seqsim(tr,length=666,model='GTR',option="normal",basefreq=rep(0.25,4))
	})

gtrGamma<-lapply(1:20,function(l){
	print(l)
seqsim(tr,length=666,model='GTR',option="discrete",basefreq=rep(0.25,4))
	})

gtrIGamma<-lapply(1:20,function(l){
	print(l)
seqsim(tr,length=666,model='GTR',option="invariant",basefreq=rep(0.25,4))->t
	saveAlignment(t,file=paste("gtrIGamma",l,".fas",sep=""),skip.internal=T)
	rm(t)
	})		


#################align_and_build_trees##################################align_and_build_trees##################################align_and_build_trees##################################align_and_build_trees##########

#align with mafft
for fasta_file in $(ls *.fas)
do
linsi $fasta_file > $fasta_file.fasta
done

#build trees
for msa_file in $(ls *.fasta)
do
./iqtree -s $msa_file
done


#################compute_spectral_density_profile##################################compute_spectral_density_profile##################################compute_spectral_density_profile################################################

treef<-list.files(path="~/Desktop/_mtree/_simulated/_nucsim/_gammatest",pattern="*.treefile")
treef<-mixedsort(treef)
trees<-lapply(treef,read.tree)


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
for(n in 1:length(trees)){
	dens(e[[n]])->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}

ped<-sapply(1:length(pe),function(l){
 	pe[l]*runif(10,0.8,1.1)
 })
skewd<-sapply(1:length(pe),function(l){
 	skew[l]*runif(10,0.8,1.1)
 })
heightd<-sapply(1:length(pe),function(l){
 	height[l]*runif(10,0.8,1.1)
 })


g<-c(rep("discrete",200),rep("invariant",200),rep("normal",200))
tab<-cbind(c(ped),c(skewd),c(heightd),g)
colnames(tab)<-c("pe","skew","height","gamma")
write.table(tab,file="Simulated_gamma_nostem_MGLtable.txt")

#################TEST_AND_PLOT##################################TEST_AND_PLOT##################################TEST_AND_PLOT#############################


##STANDARD_GRAPH_LAPLACIAN
tab<-read.table('Simulated_gamma_MGLtable.txt',header=T,row.names=1)

cols<-colors(1)[as.numeric(tab$gamma)*30]
col<-cols[seq(1,600,200)]

sactter.grid(tab$pe,tab$skew,log(tab$height),pch=20,color=cols,angle=230)

subs<-lapply(1:length(levels(tab$gamma)),function(s){
	subset(tab,tab$gamma==levels(tab$gamma)[s])
	})

pes<-list(subs[[1]]$pe, subs[[2]]$pe, subs[[3]]$pe)
skews<-list(subs[[1]]$skew, subs[[2]]$skew, subs[[3]]$skew)
heights<-list(subs[[1]]$height, subs[[2]]$height, subs[[3]]$height)

	
pdf('Simulated_gamma_phylostats_boxplot.pdf',6,12)
par(mfrow=c(3,1),mar=c(4,4,1,1))
boxplot(pes,ylab=expression(lambda),axes=F);axis(2,las=2)
stripchart(pes,vertical=TRUE,method="jitter",add=TRUE,pch=20,col= col)
boxplot(skews,ylab=expression(psi),axes=F);axis(2,las=2)
stripchart(skews,vertical=TRUE,method="jitter",add=TRUE,pch=20,col= col)
boxplot(heights,ylab=expression(eta),xlab='dN/dS',axes=F);axis(2,las=2)
	axis(1,at=1:4,label=seq(0.1,1.9,0.6))
stripchart(heights,vertical=TRUE,method="jitter",add=TRUE,pch=20,col=col) 
dev.off()

pdf('Simulated_gamma_phylostats_boxplot.pdf',6,12)
par(mfrow=c(3,2),mar=c(4,3,1,1))
color<-c(rep("brown2",20),rep("cornsilk3",20),rep("darkred",20))
lapply(seq(1,60,20),function(p){
	drop.tip(trees[[p]],1)->tree
	multi2di(tree)->trs
	plot(ladderize(trees[[p]]),show.tip.label=F,edge.color=color[p])
	plot(d[[p]]$x,dsc[[p]],type='l',xlim=c(-500,2000),ylim=c(0,0.02),
		col=color[p])
		polygon(d[[p]]$x,dsc[[p]],col=color[p])
	})
dev.off()	
	
