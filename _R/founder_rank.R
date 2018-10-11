library(igraph)
library(ape)

founder_rank_all<-function(phylos,sd=T,plot=F){
	
	#compute spectral density profile summary statistics
	gete<-function(phy){
		abs(eigen(
			graph.laplacian(
				graph.adjacency(
				data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=F),
	only.values=T)$values)
	}
	e<-lapply(phylos,gete)
	c()->d;c()->dsc;c()->pe
	for(n in 1:length(e)){
		dens(e[[n]])->d[[n]]
		d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
		max(e[[n]])->pe[[n]]
		}
	
	#define function for finding heterogeneity rank for each tree
	founder_rank<-function(phylos,test){
	#compute median+-0.5*sigma^2
	med<-function(x){
		c(median(log(x))-0.5*sd(log(x)),median(log(x))+0.5*sd(log(x)),median(log(x)))}
	
	#draw the threshold at +-0.5*sd or at the median
	if(sd==T){
		if(log(pe[[test]])<med(pe)[1]){return('homogeneous')}	
		if(log(pe[[test]])>med(pe)[2]){return('heterogeneous')}	
		else{return('medial')}
	}
	else{
		if(log(pe[[test]])<med(pe)[3]){return('homogeneous')}	
		if(log(pe[[test]])>med(pe)[3]){return('heterogeneous')}	
		}
	}
	#create table of founder ranks for all trees in set
	ranks<-sapply(1:length(phylos),function(r){
		founder_rank(phylos,r)
	})
	#return trees by PE rank and founder classification
	tab_pe<-as.data.frame(cbind(1:length(phylos),pe))
	tab<-as.data.frame(cbind(order(tab_pe[,2],tab_pe[,1]),ranks[order(tab_pe[,2],tab_pe[,1])]))
	colnames(tab)<-c('rank','founder_type')
	
	#plot PE (shifted by min) in order of tree input, demarcating +- 1/2 standard deviations from the median
	if(plot==T){
		abs(min(log(pe)))->shift
		barplot(sort((log(pe)+shift)),col=colors(1)[runif(1,30,502)],
			xlab='phylogeny',ylab=expression(lambda ~'* (shifted)'),
			names=tab$rank,cex=0.75)
		abline(h=median(log(pe))-0.5*sd(log(pe))+shift,lty=3)
		abline(h=median(log(pe))+0.5*sd(log(pe))+shift,lty=2)	
		return(tab)
	}
	else{return(tab)}
}	

