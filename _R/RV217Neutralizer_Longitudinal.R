library(RPANDA)
library(phytools)
library(igraph)
library(fpc)
library(gplots)
library(RColorBrewer)
library(seqinr)
source('~/Desktop/_mtree/_R/_JSD.R', chdir = TRUE)

trfiles<-list.files(path="~/Desktop/_mtree/RV217/AlignedByParticipant/_longitudinal",pattern="*.treefile")
	trfiles<-trfiles[-4]#remove 20263, miracle outlier
alfiles<-list.files(path="~/Desktop/_mtree/RV217/AlignedByParticipant/_longitudinal",pattern="*.fas.txt.out")
	alfiles<-alfiles[seq(1,254,9)]#reduce to only msa
	alfiles<-alfiles[-4]#remove 20263, miracle outlier
	
setwd('~/Desktop/_mtree/RV217/AlignedByParticipant/_longitudinal')
tree<-lapply(trfiles,read.tree)
align<-lapply(alfiles,read.alignment,format='fasta')

gets<-function(phy){
eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=F),
	only.values=T)$values->x
	subset(x,x>0)->e
	return(e)
}
e<-lapply(tree,gets)


c()->d;c()->dsc;c()->pe;c()->skew;c()->height
c()->gaps;c()->gapMat;c()->modalities;c()->gapMatCol;c()->eigenGap
for(n in 1:length(e)){
	dens(log(e[[n]]))->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	abs(diff(e[[n]]))->gaps[[n]]
	as.matrix(gaps[[n]])->gapMat[[n]]
	modalities[[n]] <- c(1:length(gapMat[[n]]))
    gapMatCol[[n]] <- cbind(modalities[[n]], gapMat[[n]])
    gapMatCol[[n]][order(gapMatCol[[n]][,2],decreasing=T),][1,1]->eigenGap[[n]]
	}

Ds<-c()
for(i in 1:length(e)){
	Ds<-as.data.frame(cbind(Ds,d[[i]]$x))}
	jsd<-as.matrix(JSDist(Ds-min(Ds)+1e-9))
	colnames(jsd)<-tab$bn
	rownames(jsd)<-rownames(tab)
	heatmap.2(jsd,trace="none",
			symm=T,dendrogram="column",
			cexRow=0.5,cexCol=0.5,col=brewer.pal(9,'Blues'))

	sel<-sapply(1:length(align),function(a){
	   	print(a) 
	   	kaks(align[[a]])->s
	   	if(length(s)>1){
	   	median(s$ka/s$ks,na.rm=T)
	   	}
	   	else(0)
	   	})


tab<-cbind(pe,skew,height,clus[[1]]$clustering,tips)
rownames(tab)<-files
write.csv(tabs,file="RV217neutralizer_table.csv")



#####
tab<-read.csv("RV217neutralizer_table.csv",header=T,row.names=1)

pdf('RV217longitudinal_phylospace.pdf')
sactter.grid(log(tab$pe),tab$skew,tab$height,color=tab$cluster,pch=tab$bn,angle=230,cex.symbol=2)
dev.off()	

	
subset(tab,tab$bn==1)->bn1#non-broad
subset(tab,tab$bn==2)->bn2#broad	
boxplot(log(bn1$pe),log(bn2$pe),NA,bn1$skew,bn2$skew,NA,bn1$height,bn2$height)#,ann=F,axes=F)
	axis(1,at=seq(1.5,7.5,3),labels=c("PE","skew","height"))
	axis(2,las=2)
bns<-list(log(bn1$pe),log(bn2$pe),0,bn1$skew,bn2$skew,0,bn1$height,bn2$height)
bnl<-rep(list(bn1,bn2,bn1),8)
r<-rep(c(dim(bn1)[1],dim(bn2)[1],1),8)
p<-rep(c(2,5,1),8)
s<-seq(0.75,8.75,1)
t<-seq(1.25,9.25,1)

pdf('RV217Neutralizer_broadboxplot.pdf')
#par(ann=FALSE)
plot(jitter(rep(0,10),amount=0.2), 1:10,
     xlim=range(1,8), ylim=range(-2,10),
     axes=FALSE,frame.plot=TRUE)
	lapply(1:length(bns),function(l){
		points(jitter(rep(l,r[l])),bns[[l]],pch=p[l],col=bnl[[l]]$si,cex=bnl[[l]]$mf)
		segments(s[l],mean(bns[[l]]),t[l],mean(bns[[l]]),lwd=2)
	})
	axis(1,at=seq(1.5,7.5,3),labels=c("PE","skew","height"))
	axis(2,las=2)
dev.off()



subset(tab,tab$cluster==1)->c1
subset(tab,tab$cluster==2)->c2
pdf('RV217Neutralizer_clusterboxplot.pdf')
boxplot(log(c1$pe),log(c2$pe),NA,c1$skew,c2$skew,NA,c1$height,c2$height,col=brewer.pal(9,"Blues")[c(1,9,2)],ann=F,axes=F)
	axis(1,at=seq(1.5,10.5,3),labels=c("PE","skew","height","dN/dS"))
	axis(2,las=2)
dev.off()	

pdf('RV217Neutralizer_scatterplotDimensions.pdf',4,10)
par(mfrow=c(2,1))
plot(tab$kaks~tab$skew,col=c("#A6CEE3","#08306B")[tab$cluster],pch=tab$bn,xlim=c(0,11),cex=2)
  abline(lm(tab$kaks~tab$skew))
  legend('topleft',c('y=1+0.06x,r2=0.25,p=0.003'),bty="n",lty=1)
plot(log(tab$tips)~log(tab$pe),col=c("#A6CEE3","#08306B")[tab$cluster],pch=tab$bn,cex=2)
  abline(lm(log(tab$tips)~log(tab$pe)))
  legend('topleft',c('y=3.79+0.27x,r2=0.52,p<0.001'),bty="n",lty=1)
dev.off()  
  
c()->t
for(s in 1:length(levels(tab$subtype))){
	subset(tab,tab$subtype==levels(tab$subtype)[s])->t[[s]]
}

boxplot(log(t[[1]]$pe),log(t[[2]]$pe),log(t[[3]]$pe),log(t[[4]]$pe))
boxplot(t[[1]]$skew,t[[2]]$skew,t[[3]]$skew,t[[4]]$skew)
boxplot(t[[1]]$height,t[[2]]$height,t[[3]]$height,t[[4]]$height)
