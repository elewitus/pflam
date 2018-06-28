library(RPANDA)
library(phytools)
library(igraph)
library(gplots)
library(RColorBrewer)
library(fpc)
source('~/Desktop/_mtree/_R/_JSD.R')

setwd('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/_genemsa')
folders<-list.files(path="~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/_genemsa",pattern="")
files<-lapply(1:length(folders),function(f){
	list.files(path=paste("~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/_genemsa/",folders[f],sep=""),pattern="*.treefile")
	})	
c()->trees	
trees<-lapply(1:length(files),function(t){
	print(t)
	setwd(paste("~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/_genemsa/",folders[t],sep=""))
	lapply(files[[t]],read.tree)
	})

filenames<-c(files[[1]],files[[2]],files[[3]],files[[4]],files[[5]],files[[6]],files[[7]],files[[8]],files[[9]])
genetrees<-c(trees[[1]],trees[[2]],trees[[3]],trees[[4]],trees[[5]],trees[[6]],trees[[7]],trees[[8]],trees[[9]])
###########standard_spectral_profiles#####################standard_spectral_profiles#####################standard_spectral_profiles#####################standard_spectral_profiles#####################
c()->all_eigen
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
	rep(paste(folders[r]),length(trees[[r]]))
		})
		c(rs[,1],rs[,2],rs[,3],rs[,4],rs[,5],rs[,6],rs[,7],rs[,8],rs[,9])->rows
		c(rep("S",172),rep("R",86),rep("A",129))->type
		
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
	colnames(virus_JSD)<-type
heatmap.2(virus_JSD,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.3,cexCol=0.3,col=cols)

pamk.best <- pamk(virus_JSD)

#multidimensional plot
c()->pe;c()->height;c()->skew;c()->tips
for(q in 1:length(d)){
	max(x[[q]])->pe[[q]]
	log(max(dsc[[q]]))->height[[q]]
	log(skewness(dsc[[q]]))->skew[[q]]
	Ntip(genetrees[[q]])->tips[[q]]
	}
	tab<-cbind(pe,height,skew,type,rows,filenames,pamk.best[[1]]$clustering,tips)
		colnames(tab)<-c("pe","height","skew","type","gene","participant","cluster","tips")
		rownames(tab)<-1:dim(tab)[1]
	write.table(tab,file="RV217_AlignedByParticipant_genes_MGLtable.txt")
	
###########plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################plotting_spectral_density_profile#####################	##########plotting_spectral_density_profile##############################
colorbyfunction<-brewer.pal(3,"Paired")
cbf<-brewer.pal(3,"Paired")[c(1,1,1,1,2,2,3,3,3)]
colorbygene<-brewer.pal(9,"Paired")

tab<-read.table('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/RV217_AlignedByParticipant_genes_MGLtable.txt',header=T,row.names=1)

pdf('RV217_AlignedByParticipant_genes_phylospace.pdf')
sactter.grid(tab$pe,tab$height,tab$skew,angle=30,xlab="PE",ylab="height",zlab="skew",pch=20,color=colorbyfunction)
dev.off()

tabbygene<-lapply(1:length(folders),function(tbg){
	subset(tab,tab$gene==paste(folders[tbg]))
	})
	
pdf('RV217_AlignedByParticipant_genes_phylostats.pdf',10,8)
par(mfcol=c(3,2),mar=c(4,4,1,1))
	boxplot(tabbygene[[1]]$pe,tabbygene[[2]]$pe,tabbygene[[3]]$pe,tabbygene[[4]]$pe,tabbygene[[5]]$pe,tabbygene[[6]]$pe,tabbygene[[7]]$pe,tabbygene[[8]]$pe,tabbygene[[9]]$pe,outline=F,col= cbf,ann=F,ylab="PE")
	boxplot(tabbygene[[1]]$skew,tabbygene[[2]]$skew,tabbygene[[3]]$skew,tabbygene[[4]]$skew,tabbygene[[5]]$skew,tabbygene[[6]]$skew,tabbygene[[7]]$skew,tabbygene[[8]]$skew,tabbygene[[9]]$skew,outline=F,col= cbf,ann=F,ylab="skewness")
	boxplot(tabbygene[[1]]$height,tabbygene[[2]]$height,tabbygene[[3]]$height,tabbygene[[4]]$height,tabbygene[[5]]$height,tabbygene[[6]]$height,tabbygene[[7]]$height,tabbygene[[8]]$height,tabbygene[[9]]$height,outline=F,col= cbf,names=folders,las=2,ylab="height")

for(o in 1:length(tabbygene)){
	print(o)
	plot(density(tabbygene[[o]]$pe),col= cbf[o],ann=F,axes=F,xlim=c(-1,5),ylim=c(0,3),lwd=2)
	par(new=T)
	}
	axis(1)
	axis(2)
	par(new=F)
	legend("topright",c("structural","regulatory","accessory"),col=brewer.pal(3,"Paired"),lty=1,bty="n")
for(o in 1:length(tabbygene)){
	plot(density(tabbygene[[o]]$skew),col= cbf[o],axes=F,ann=F,xlim=c(0,7),ylim=c(0,2),lwd=2)
	par(new=T)
	}
	axis(1)
	axis(2)
	par(new=F)
for(o in 1:length(tabbygene)){
	plot(density(tabbygene[[o]]$height),col=cbf[o],axes=F,ann=F,xlim=c(-5,18),ylim=c(0,0.5),lwd=2)
	par(new=T)
	}
	axis(1)
	axis(2)		
dev.off()

###############skew_effect#####################skew_effect##############################skew_effect##############################skew_effect####################################################################


bigskew<-lapply(5:length(tabbygene),function(r){
	subset(tabbygene[[r]],tabbygene[[r]]$skew>3.5)$participant
	})
smallskew<-lapply(1:4,function(r){
	subset(tabbygene[[r]],tabbygene[[r]]$skew<3)$participant
	})	

c()->overlap
overlaps<-sapply(1:length(bigskew),function(a){
overlap[[a]]<-sapply(1:length(bigskew),function(b){
	length(intersect(bigskew[[a]],bigskew[[b]]))/length(bigskew[[b]])
	})
})

overlaps<-lapply(1:length(bigskew),function(a){
overlap<-lapply(1:length(bigskew),function(b){
	intersect(bigskew[[a]],bigskew[[b]])})
})

Reduce(intersect, list(bigskew[[1]],bigskew[[2]],bigskew[[3]],bigskew[[4]],bigskew[[5]]))->RAgenes
Reduce(intersect, list(smallskew[[1]],smallskew[[2]],smallskew[[3]],smallskew[[4]]))->Sgenes
intersect(RAgenes,Sgenes)

dif<-lapply(1:4,function(d){
	print(d)
	setdiff(vargene,bigskew[[d]])->vargene
})

Reduce(intersect, list(dif[[1]], dif[[3]], dif[[4]]))

corrgram(overlaps,labels=folders,lower.panel=panel.shade,upper.panel=NULL)

par(mfrow=c(3,3))
lapply(1:length(trees),function(z){
	tabs<-subset(tab,tab$gene==paste(folders[z]))
	#tabs<-tab[c((seq(1,387,42)[z]+1):(seq(1,387,42)[z+1]+1)),]
sactter.grid(tabs$pe,log(tabs$height),log(tabs$skew),angle=30,xlab="PE",ylab="height",zlab="skew",color=as.numeric(tabs$gene),pch=20,main=paste(folders[z]))
})


###############multifounder_effect#####################multifounder_effect##############################multifounder_effect##############################multifounder_effect####################################
colorbyfunction<-brewer.pal(3,"Paired")
cbf<-brewer.pal(3,"Paired")[c(1,1,1,1,2,2,3,3,3)]
colorbygene<-brewer.pal(9,"Paired")


tab<-read.table('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/RV217_AlignedByParticipant_genes_MGLtable_noLog.txt',header=T,row.names=1)

multif<-read.table('~/Desktop/_mtree/RV217/multifounder.txt',header=T)
	multif<-multif[-7,]
mfs<-which(tab$participant%in%multif)

##multifounder
tab[mfs,]->mftab

sactter.grid(mftab$pe,mftab$height,mftab$skew,angle=30,xlab="PE",ylab="height",zlab="skew",pch=20,color=colorbyfunction[as.numeric(mftab$type)])

mfbygene<-lapply(1:length(folders),function(tbg){
	subset(mftab,mftab$gene==paste(folders[tbg]))
	})

par(mfcol=c(3,2),mar=c(4,4,1,1))
	boxplot(mfbygene[[1]]$pe,mfbygene[[2]]$pe,mfbygene[[3]]$pe,mfbygene[[4]]$pe,mfbygene[[5]]$pe,mfbygene[[6]]$pe,mfbygene[[7]]$pe,mfbygene[[8]]$pe,mfbygene[[9]]$pe,outline=F,col= cbf,ann=F,ylab="PE")
	boxplot(mfbygene[[1]]$skew,mfbygene[[2]]$skew,mfbygene[[3]]$skew,mfbygene[[4]]$skew,mfbygene[[5]]$skew,mfbygene[[6]]$skew,mfbygene[[7]]$skew,mfbygene[[8]]$skew,mfbygene[[9]]$skew,outline=F,col= cbf,ann=F,ylab="skewness")
	boxplot(mfbygene[[1]]$height,mfbygene[[2]]$height,mfbygene[[3]]$height,mfbygene[[4]]$height,mfbygene[[5]]$height,mfbygene[[6]]$height,mfbygene[[7]]$height,mfbygene[[8]]$height,mfbygene[[9]]$height,outline=F,col= cbf,names=folders,las=2,ylab="height")

#cohort_effect#female cohorts - UG(10XXX), KY(20XXX), TZ(30XXX)
#male cohort - TH(40XXX)
colorbycohort<-c(brewer.pal(3,"Blues"),"#E41A1C")
mfc<-lapply(1:4,function(e){
	subset(mftab,mftab$cohort==levels(mftab$cohort)[e])
	})
	pdf('RV217_AlignedByParticipant_genes_founder_gender_and_cohort.pdf')
	par(mfcol=c(3,2),mar=c(2,2,1,1))
	boxplot(mfc[[1]]$pe,mfc[[4]]$pe,mfc[[3]]$pe,mfc[[2]]$pe,col=colorbycohort,outline=F, ylab ="PE")
	boxplot(mfc[[1]]$skew,mfc[[4]]$skew,mfc[[3]]$skew,mfc[[2]]$skew,col=colorbycohort,outline=F,ylab="skew")
	boxplot(mfc[[1]]$height,mfc[[4]]$height,mfc[[3]]$height,mfc[[2]]$height,col= colorbycohort,outline=F, ylab ="height",names=c("KY","UG","TZ","TH"))
		


#par(mfcol=c(3,9),mar=c(2,2,1,1))
#lapply(1:length(mfbygene),function(c){
#	barplot(mfbygene[[c]]$pe,col=colorbycohort[as.numeric(mfbygene[[c]]$cohort)])
#	barplot(mfbygene[[c]]$skew,col=colorbycohort[as.numeric(mfbygene[[c]]$cohort)])
#	barplot(mfbygene[[c]]$height,col=colorbycohort[as.numeric(mfbygene[[c]]$cohort)])		
#	})

##singlefounder
tab[-mfs,]->sftab

sactter.grid(sftab$pe,sftab$height,sftab$skew,angle=30,xlab="PE",ylab="height",zlab="skew",pch=20,color=colorbyfunction[as.numeric(sftab$type)])

sfbygene<-lapply(1:length(folders),function(tbg){
	subset(sftab,sftab$gene==paste(folders[tbg]))
	})

	boxplot(sfbygene[[1]]$pe,sfbygene[[2]]$pe,sfbygene[[3]]$pe,sfbygene[[4]]$pe,sfbygene[[5]]$pe,sfbygene[[6]]$pe,sfbygene[[7]]$pe,sfbygene[[8]]$pe,sfbygene[[9]]$pe,outline=F,col= cbf,ann=F,ylab="")
	boxplot(sfbygene[[1]]$skew,sfbygene[[2]]$skew,sfbygene[[3]]$skew,sfbygene[[4]]$skew,sfbygene[[5]]$skew,sfbygene[[6]]$skew,sfbygene[[7]]$skew,sfbygene[[8]]$skew,sfbygene[[9]]$skew,outline=F,col= cbf,ann=F,ylab="")
	boxplot(sfbygene[[1]]$height,sfbygene[[2]]$height,sfbygene[[3]]$height,sfbygene[[4]]$height,sfbygene[[5]]$height,sfbygene[[6]]$height,sfbygene[[7]]$height,sfbygene[[8]]$height,sfbygene[[9]]$height,outline=F,col= cbf,names=folders,las=2,ylab="")

#cohort_effect
sfc<-lapply(1:4,function(e){
	subset(sftab,sftab$cohort==levels(sftab$cohort)[e])
	})
	
	boxplot(sfc[[1]]$pe,sfc[[4]]$pe,sfc[[3]]$pe,sfc[[2]]$pe,col=colorbycohort,outline=F)
	boxplot(sfc[[1]]$skew,sfc[[4]]$skew,sfc[[3]]$skew,sfc[[2]]$skew,col=colorbycohort,outline=F)
	boxplot(sfc[[1]]$height,sfc[[4]]$height,sfc[[3]]$height,sfc[[2]]$height,col= colorbycohort,outline=F,names=c("KY","UG","TZ","TH"))
	dev.off()

	
###############SPVL###########SPVL##########SPVL##########SPVL##############SPVL#########SPVL############SPVL###################SPVL####################SPVL#####SPVL##############SPVL##########SPVL##########

tab<-read.csv('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/RV217_AlignedByParticipant_genes_MGLtable.csv',header=T,row.names=1)


mfpe<-sapply(1:length(mfbygene),function(s){
summary(lm(mfbygene[[s]]$spvl~ mfbygene[[s]]$pe))[[9]]})
mfskew<-sapply(1:length(mfbygene),function(s){
summary(lm(mfbygene[[s]]$spvl~ mfbygene[[s]]$skew))[[9]]})
mfheight<-sapply(1:length(mfbygene),function(s){
summary(lm(mfbygene[[s]]$spvl~ mfbygene[[s]]$height))[[9]]})
cbind(mfpe,mfskew,mfheight)


sfpe<-sapply(1:length(sfbygene),function(s){
summary(lm(sfbygene[[s]]$spvl~ sfbygene[[s]]$pe))[[9]]})
sfskew<-sapply(1:length(sfbygene),function(s){
summary(lm(sfbygene[[s]]$spvl~ sfbygene[[s]]$skew))[[9]]})
sfheight<-sapply(1:length(sfbygene),function(s){
summary(lm(sfbygene[[s]]$spvl~ sfbygene[[s]]$height))[[9]]})
cbind(sfpe,sfskew,sfheight)

