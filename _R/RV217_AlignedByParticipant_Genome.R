library(RPANDA)
library(phytools)
library(igraph)
library(gplots)
library(RColorBrewer)
library(fpc)
source('~/Desktop/_mtree/_R/_JSD.R')
cols<-brewer.pal(9,"Reds")
col <- colorRampPalette(cols)(n = 20)

setwd('~/Desktop/_mtree/RV217/AlignedByParticipant/Genome')
list.files(path="~/Desktop/_mtree/RV217/AlignedByParticipant/Genome",pattern="*.treefile")->files
files<-files[-c(7,43)]
trees<-lapply(files,read.tree)

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
        subset(x[[i]], x[[i]] >= 0) -> x[[i]]  		
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
	rownames(virus_JSD)<-files
	colnames(virus_JSD)<-files
heatmap.2(virus_JSD,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.5,cexCol=0.5,col=cols)

pamk.best <- pamk(virus_JSD)

#multidimensional plot
c()->pe;c()->height;c()->skew;c()->tips
for(q in 1:length(d)){
	log(max(x[[q]]))->pe[[q]]
	max(dsc[[q]])->height[[q]]
	skewness(dsc[[q]])->skew[[q]]
	}
	tab<-cbind(pe,height,skew,pamk.best[[1]]$clustering)
		colnames(tab)<-c("pe","height","skew","cluster")

	
###########plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################	##########plotting_spectral_density_profile##############################

tab<-read.table('~/Desktop/_mtree/RV217/AlignedByParticipant/Genome/RV217_AlignedByParticipant_Genome_MGLtable.txt',header=T,row.names=1)

pdf('RV217_AlignedByParticipant_phylospace.pdf',10,5)
par(mfrow=c(1,3))
plot(pam(virus_JSD,pamk.best$nc))
sactter.grid(tab$pe,tab$height,tab$skew,angle=30,xlab="PE",ylab="height",zlab="skew",color=1,pch=as.numeric(tab$multifounder))
dev.off()			

