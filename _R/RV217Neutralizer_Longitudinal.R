library(RPANDA)
library(phytools)
library(igraph)
library(fpc)
library(gplots)
library(RColorBrewer)
source('~/Desktop/_mtree/_R/_JSD.R', chdir = TRUE)

files<-list.files(path="~/Desktop/_mtree/RV217/AlignedByParticipant/_longitudinal",pattern="*.treefile")
	files<-files[-4]
setwd('~/Desktop/_mtree/RV217/AlignedByParticipant/_longitudinal')
tree<-lapply(files,read.tree)

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
for(n in 1:length(e)){
	dens(log(e[[n]]))->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
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

	
tab<-cbind(pe,skew,height,clus[[1]]$clustering,tips)
rownames(tab)<-files
write.csv(tabs,file="RV217neutralizer_table.csv")



#####
tab<-read.csv("RV217neutralizer_table.csv",header=T,row.names=1)

pdf('RV217longitudinal_phylospace.pdf')
sactter.grid(log(tab$pe),tab$skew,tab$height,color=brewer.pal(3,"Paired")[tab$cluster],pch=tab$bn+16,angle=230,cex.symbol=2)
dev.off()	
	
	
subset(tab,tab$bn==1)->bn1#non-broad
subset(tab,tab$bn==2)->bn2#broad	
boxplot(log(bn1$pe),log(bn2$pe),NA,bn1$skew,bn2$skew,NA,bn1$height,bn2$height)#,ann=F,axes=F)
	axis(1,at=seq(1.5,7.5,3),labels=c("PE","skew","height"))
	axis(2,las=2)
bns<-list(log(bn1$pe),log(bn2$pe),0,bn1$skew,bn2$skew,0,bn1$height,bn2$height)
bnl<-rep(list(bn1,bn2,bn1),8)
r<-rep(c(9,15,1),8)
p<-rep(c(2,5,1),8)
s<-seq(0.75,8.75,1)
t<-seq(1.25,9.25,1)

pdf('RV217Neutralizer_broadboxplot.pdf')
#par(ann=FALSE)
plot(jitter(rep(0,10),amount=0.2), t1,
     xlim=range(1,8), ylim=range(-2,10),
     axes=FALSE,frame.plot=TRUE)
	lapply(1:length(bns),function(l){
		points(jitter(rep(l,r[l])),bns[[l]],pch=p[l],col=bnl[[l]]$si,cex=bnl[[l]]$mf)
		segments(s[l],mean(bns[[l]]),t[l],mean(bns[[l]]),lwd=2)
	})
	axis(1,at=seq(1.5,7.5,3),labels=c("PE","skew","height"))
	axis(2,las=2)
dev.off()



subset(tabs,tabs$cluster==1)->c1
subset(tabs,tabs$cluster==2)->c2
pdf('RV217Neutralizer_clusterboxplot.pdf')
boxplot(log(c1$pe),log(c2$pe),NA,c1$skew,c2$skew,NA,c1$height,c2$height,col=cols[c(1,9,2)],ann=F,axes=F)
	axis(1,at=seq(1.5,7.5,3),labels=c("PE","skew","height"))
	axis(2,las=2)
dev.off()	
