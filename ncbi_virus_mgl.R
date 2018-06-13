library(RPANDA)
library(phytools)
library(igraph)
library(gplots)
library(RColorBrewer)
source('_JSD.R')
cols<-brewer.pal(9,"Blues")
col <- colorRampPalette(cols)(n = 20)

setwd('~/Desktop/_newClades/_virus/_tr')
files <- list.files(path="~/Desktop/_newClades/_virus/_tr",pattern="*.treefile")
trees <- lapply(files, read.tree)
names<-read.table('virus_tree_names.txt',header=T,row.names=1)


###########standard_spectral_profiles#####################standard_spectral_profiles#####################standard_spectral_profiles#####################standard_spectral_profiles#####################
	trAll <- lapply(trees,dist.nodes)	 
	allMats <- lapply(trAll,data.matrix)
	all_graph <- lapply(allMats,graph.adjacency,weighted=T)
	all_laplacian <- lapply(all_graph,graph.laplacian,
					normalized=F)			
	all_eigen <- lapply(all_laplacian,eigen,
				symmetric=TRUE,only.values=TRUE)	
				
##prepare for convolution
list() -> x; list() -> n
for(i in 1:length(all_eigen))
        {
        abs(all_eigen[[i]]$values) -> x[[i]]
        subset(x[[i]], x[[i]] >= 1) -> x[[i]]  		
   		length(x[[i]]) -> n[[i]]
        }

##Convolve with gaussian kernel, f(x)/integrate{f(y)}
list() -> d; list() -> dint; list() -> dsc
for(i in 1:length(x)){
	dens(log(x[[i]])) -> d[[i]]
	integr(d[[i]]$x,d[[i]]$y) -> dint[[i]]
	(d[[i]]$y/dint[[i]]) -> dsc[[i]]
	}  

#clustering
Ds<-c() 
for(i in 1:length(x)){
     Ds<-as.data.frame(cbind(Ds,abs(d[[i]]$x)))
     }
virus_JSD <- as.matrix(JSDist(Ds))
	rownames(virus_JSD)<-names[,1]
	colnames(virus_JSD)<-names[,1]
heatmap.2(virus_JSD,,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.5,cexCol=0.5,col=cols)

pamk.best <- pamk(virus_JSD)

#multidimensional plot
c()->pe;c()->height;c()->skew;c()->tips
for(q in 1:length(d)){
	log(max(x[[q]]))->pe[[q]]
	max(dsc[[q]])->height[[q]]
	skewness(dsc[[q]])->skew[[q]]
	log(Ntip(trees[[q]]))->tips[[q]]
	}
	tab<-cbind(pe,height,skew,tips,pamk.best[[1]]$clustering)
		colnames(tab)<-c("pe","height","skew","tips","cluster")

sactter.grid(pe,height,skew,pch=21,bg="white",angle=230,xlab="PE",ylab="height",zlab="skew",color=pamk.best[[1]]$clustering)

	c1<-subset(tab,tab[,5]==1)
	c2<-subset(tab,tab[,5]==2)	
boxplot(c1[,1],c2[,1],NA,c1[,2],c2[,2],NA,c1[,3],c2[,3],NA,c1[,4],c2[,4],ann=F)
	axis(1,at=seq(1.5,10.5,3),labels=colnames(tab)[1:4])
	axis(2,las=2)

##################normal_spectral_profiles##################normal_spectral_profiles##################normal_spectral_profiles##################normal_spectral_profiles###############################
	trAll <- lapply(trees,dist.nodes)	 
	allMats <- lapply(trAll,data.matrix)
	all_graph <- lapply(allMats,graph.adjacency,weighted=T)
	all_laplacian <- lapply(all_graph,graph.laplacian,
					normalized=T)			
	all_eigen <- lapply(all_laplacian,eigen,
				symmetric=TRUE,only.values=TRUE)	
				
##prepare for convolution
list() -> x; list() -> n
for(i in 1:length(all_eigen))
        {
        abs(all_eigen[[i]]$values) -> x[[i]]
        subset(x[[i]], x[[i]] >= 1) -> x[[i]]  		
   		length(x[[i]]) -> n[[i]]
        }

##Convolve with gaussian kernel, f(x)/integrate{f(y)}
list() -> d; list() -> dint; list() -> dsc
for(i in 1:length(x)){
	dens(log(x[[i]])) -> d[[i]]
	integr(d[[i]]$x,d[[i]]$y) -> dint[[i]]
	(d[[i]]$y/dint[[i]]) -> dsc[[i]]
	}  

#clustering
Ds<-c() 
for(i in 1:length(x)){
     Ds<-as.data.frame(cbind(Ds,abs(d[[i]]$x)))
     }
virus_JSD <- as.matrix(JSDist(Ds))
	rownames(virus_JSD)<-names[,1]
	colnames(virus_JSD)<-tab[,5]
heatmap.2(virus_JSD,,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.5,cexCol=0.5,col=cols)

pamk.best <- pamk(virus_JSD)

#multidimensional plot
c()->pe;c()->height;c()->skew;c()->tips
for(q in 1:length(d)){
	max(x[[q]])->pe[[q]]
	log(max(dsc[[q]]))->height[[q]]
	log(skewness(dsc[[q]]))->skew[[q]]
	log(Ntip(trees[[q]]))->tips[[q]]
	}
	tab<-cbind(pe,height,skew,tips,pamk.best[[1]]$clustering)
		colnames(tab)<-c("pe","height","skew","tips","cluster")

sactter.grid(pe,height,skew,pch=21,bg="white",angle=230,xlab="PE",ylab="height",zlab="skew",color=pamk.best[[1]]$clustering)

	c1<-subset(tab,tab[,5]==1)
	c2<-subset(tab,tab[,5]==2)
	c3<-subset(tab,tab[,5]==3)	
boxplot(c1[,1],c2[,1],c3[,1],NA,c1[,2],c2[,2],c3[,2],NA,c1[,3],c2[,3],c3[,3],NA,c1[,4],c2[,4],c3[,4],ann=F,axes=F)
	axis(1,at=seq(2,16,4),labels=colnames(tab)[1:4])
	axis(2,las=2)

plot(d[[1]]$x,dsc[[1]],type="l",col=tab[1,5])
sapply(1:length(d),function(q){
	plot(d[[q]]$x,log(dsc[[q]]),
	type="l",col=tab[q,5],
	xlim=c(0,1),ylim=c(-1e3,0),
	ann=F,axes=F)
	par(new=T)
})
