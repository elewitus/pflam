library(RPANDA)
library(phytools)
library(igraph)
library(gplots)
library(RColorBrewer)
library(fpc)
source('~/Desktop/_mtree/_R/_JSD.R')
cols<-brewer.pal(9,"Reds")
col <- colorRampPalette(cols)(n = 20)

setwd('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes')
list.files(path="~/Desktop/_mtree/RV217/AlignedByParticipant/_genes",pattern="*")->folders
files<-lapply(1:length(folders),function(f){
	list.files(path=paste("~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/",folders[f],sep=""),pattern="*.treefile")
	})
c()->trees	
trees<-lapply(1:length(files),function(t){
	setwd(paste("~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/",folders[t],sep=""))
	lapply(files[[t]],read.tree)
	})
#genetrees<-c(trees[[1]],trees[[2]],trees[[3]],trees[[4]],trees[[5]],trees[[6]],trees[[7]],trees[[8]],trees[[9]])
###########standard_spectral_profiles#####################standard_spectral_profiles#####################standard_spectral_profiles#####################standard_spectral_profiles#####################

for(u in 1:length(trees)){	
	trAll <- lapply(trees[[u]],dist.nodes)	 
	allMats <- lapply(trAll,data.matrix)
	all_graph <- lapply(allMats,graph.adjacency,weighted=T)
	all_laplacian <- lapply(all_graph,graph.laplacian,
					normalized=F)			
	all_eigen[[u]] <- lapply(all_laplacian,eigen,
				symmetric=TRUE,only.values=TRUE)	
}		
all_eigens<-c(all_eigen[[1]],all_eigen[[2]],all_eigen[[3]],all_eigen[[4]],all_eigen[[5]],all_eigen[[6]],all_eigen[[7]],all_eigen[[8]],all_eigen[[9]])

rs<-sapply(1:length(folders),function(r){
	rep(paste(folders[r]),43)
		})
		c(rs[,1],rs[,2],rs[,3],rs[,4],rs[,5],rs[,6],rs[,7],rs[,8],rs[,9])->rows
		
##prepare for convolution
list() -> x; list() -> n
for(i in 1:length(all_eigens))
        {
        abs(all_eigens[[i]]$values) -> x[[i]]
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
	rownames(virus_JSD)<-rows
	colnames(virus_JSD)<-rows
heatmap.2(virus_JSD,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.3,cexCol=0.3,col=cols)

pamk.best <- pamk(virus_JSD,krange=3)

#multidimensional plot
c()->pe;c()->height;c()->skew;c()->tips
for(q in 1:length(d)){
	max(x[[q]])->pe[[q]]
	log(max(dsc[[q]]))->height[[q]]
	log(skewness(dsc[[q]]))->skew[[q]]
	}
	tab<-cbind(pe,height,skew,rows, filenames,pamk.best[[1]]$clustering)
		colnames(tab)<-c("pe","height","skew","gene","participant","cluster")
		rownames(tab)<-1:dim(tab)[1]
	write.table(tab,file="RV217_AlignedByParticipant_genes_MGLtable.txt")
	
###########plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################	##########plotting_spectral_density_profile##############################
colors<-brewer.pal(9,"Paired")
tab<-read.table('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/RV217_AlignedByParticipant_genes_MGLtable.txt',header=T,row.names=1)

pdf('RV217_AlignedByParticipant_phylospace.pdf')
#par(mfrow=c(1,2))
#plot(pam(virus_JSD,pamk.best$nc))
sactter.grid(tab$pe,tab$height,tab$skew,angle=30,xlab="PE",ylab="height",zlab="skew",color=colors[as.numeric(tab$gene)],pch=20)
dev.off()			

tabbygene<-lapply(1:length(trees),function(tbg){
	subset(tab,tab$gene==paste(folders[tbg]))
	})

par(mfcol=c(3,2),mar=c(4,4,1,1))
	boxplot(tabbygene[[1]]$pe,tabbygene[[2]]$pe,tabbygene[[3]]$pe,tabbygene[[4]]$pe,tabbygene[[5]]$pe,tabbygene[[6]]$pe,tabbygene[[7]]$pe,tabbygene[[8]]$pe,tabbygene[[9]]$pe,outline=F,col=colors,ann=F,ylab="PE")
	boxplot(tabbygene[[1]]$skew,tabbygene[[2]]$skew,tabbygene[[3]]$skew,tabbygene[[4]]$skew,tabbygene[[5]]$skew,tabbygene[[6]]$skew,tabbygene[[7]]$skew,tabbygene[[8]]$skew,tabbygene[[9]]$skew,outline=F,col=colors,ann=F,ylab="skewness")
	boxplot(tabbygene[[1]]$height,tabbygene[[2]]$height,tabbygene[[3]]$height,tabbygene[[4]]$height,tabbygene[[5]]$height,tabbygene[[6]]$height,tabbygene[[7]]$height,tabbygene[[8]]$height,tabbygene[[9]]$height,outline=F,col=colors,names=folders,las=2,ylab="height")

for(o in 1:length(files)){
	plot(density(tabbygene[[o]]$pe),col=colors[o],ann=F,axes=F,xlim=c(-1,5),ylim=c(0,3),lwd=2)
	par(new=T)
	}
	axis(1)
	axis(2)
	par(new=F)
for(o in 1:length(files)){
	plot(density(tabbygene[[o]]$skew),col=colors[o],axes=F,ann=F,xlim=c(0,7),ylim=c(0,2),lwd=2)
	par(new=T)
	}
	axis(1)
	axis(2)
	par(new=F)
for(o in 1:length(files)){
	plot(density(tabbygene[[o]]$height),col=colors[o],axes=F,ann=F,xlim=c(-2,10),ylim=c(0,2),lwd=2)
	par(new=T)
	}
	axis(1)
	axis(2)		


bigskew<-lapply(1:length(tabbygene),function(r){
	subset(tabbygene[[r]],tabbygene[[r]]$skew>3.5)$participant
	})

c()->overlap
overlaps<-sapply(1:length(bigskew),function(a){
overlap[[a]]<-sapply(1:length(bigskew),function(b){
	length(which(bigskew[[a]]%in%bigskew[[b]]))/length(bigskew[[a]])
	})
	})

corrgram(overlaps,labels=folders,lower.panel=panel.shade,upper.panel=NULL)

par(mfrow=c(3,3))
lapply(1:length(trees),function(z){
	tabs<-subset(tab,tab$gene==paste(folders[z]))
	#tabs<-tab[c((seq(1,387,42)[z]+1):(seq(1,387,42)[z+1]+1)),]
sactter.grid(tabs$pe,log(tabs$height),log(tabs$skew),angle=30,xlab="PE",ylab="height",zlab="skew",color=as.numeric(tabs$gene),pch=20,main=paste(folders[z]))
})
