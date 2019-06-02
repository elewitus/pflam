library(igraph)
library(ape)
library(fpc)

founder_rank<-function(phylos,sd=F,test=c('median','jump','partition'),g=0,plot=F){
	
#gaussian kernel
sigma = 0.1
gKernel <- function(x) 1/(sigma*sqrt(2*pi)) * exp(-(x^2)/2*sigma^2)
kernelG <- function(x, mean=0, sd=1) dnorm(x, mean = mean, sd = sd)

#kernel density estimate
dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    x <- na.omit(x)
    if(length(x) == 0)
        stop("no finite or non-missing data!")
  }
  sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

integr <- function(x, f)
{
       
       # var is numeric
       if (!is.numeric(x))
       {
              stop('The variable of integration "x" is not numeric.')
       }

       # integrand is numeric
       if (!is.numeric(f))
       {
              stop('The integrand "f" is not numeric.')
       }

       # length(var)=length(int)
       if (length(x) != length(f))
       {
              stop('The lengths of the variable of integration and the integrand do not match.')
       }

      # get lengths of var and integrand
       n = length(x)

       # trapezoidal integration
       integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))

       # print definite integral
       return(integral)
}
	
	#compute graph Laplacian eigenvalues
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
	
	#compute profile, get principle eigenvalue
	c()->d;c()->dsc;c()->pe
	for(n in 1:length(e)){
		dens(e[[n]])->d[[n]]
		d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
		max(e[[n]])->pe[[n]]
		}
		
	#define function for finding heterogeneity rank for each tree
	median_founder_test<-function(phylos,test){
	
	#compute median+-0.5*sigma^2
	g=g
	med<-function(x){
		c(median(log(x))-g*sd(log(x)),median(log(x))+g*sd(log(x)),median(log(x)))}
	
	#draw the threshold at +-g*sd or at the median
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
	
	jump_founder_test<-function(pe){
		sortedpe<-sort(log(pe))[-c(1,length(pe))]
		which.max(diff(sortedpe))->threshold
		sapply(1:length(pe),function(o){
			if(log(pe[[o]])<pe[[threshold]]){return('homogeneous')}
			if(log(pe[[o]])>pe[[threshold]]){return('heterogeneous')}
		})
	}
	
	clus_founder_test<-function(pe){
		pamk(pe,2)->clus
		sapply(1:length(clus[[1]]$partitioning),function(h){
			if(clus[[1]]$partitioning[h]==1){return('homogeneous')}
			if(clus[[1]]$partitioning[h]==2){return('heterogeneous')}
		})
	}
	
	#create table of founder ranks for all trees in set
	if(test=='median'){
	ranks<-sapply(1:length(phylos),function(r){
		median_founder_test(phylos,r)
	})
	}
	
	if(test=='jump'){
		ranks<-jump_founder_test(pe)
	}
	
	if(test=='partition'){
		ranks<-clus_founder_test(pe)
	}
	
	#return trees by PE rank and founder classification
	tab_pe<-as.data.frame(cbind(1:length(phylos),pe))
	tab<-as.data.frame(cbind(as.numeric(order(tab_pe[,2],tab_pe[,1])),ranks[order(tab_pe[,2],tab_pe[,1])]))
	colnames(tab)<-c('tree','founder_type')
	
	#plot PE (shifted by min) in order of tree input, demarcating +- 1/2 standard deviations from the median
	if(plot==T & test=='median'){
		min(log(pe))->shift
		barplot(sort((log(pe)-shift)),col=colors(1)[runif(1,30,502)],
			xlab='phylogeny',ylab=expression('ln'~lambda ~'* (shifted)'),
			names=tab$rank,cex.names=0.75)
		abline(h=median(log(pe))-g*sd(log(pe))-shift,lty=3)
		abline(h=median(log(pe))+g*sd(log(pe))-shift,lty=2)	
		abline(h=median(log(pe)),lty=1)
		return(tab)
	}
	if(plot==T & test=='jump'){
		min(log(pe))->shift
		barplot(sort((log(pe)-shift)),col=colors(1)[runif(1,30,502)],
			xlab='phylogeny',ylab=expression('ln'~lambda ~'* (shifted)'),
			names=tab$rank,cex.names=0.75)
			sortedpe<-sort(log(pe))[-c(1,length(pe))]
			sortedpe[which.max(diff(sortedpe))]->threshold
		abline(h=threshold-shift,lty=3)
		return(tab)
	}
	if(plot==T & test=='partition'){
		min(log(pe))->shift
		barplot(sort(log(pe)-shift),col=colors(1)[runif(1,30,502)],
			xlab='phylogeny',ylab=expression('ln'~lambda ~'* (shifted)'),
			names=tab$rank,cex.names=0.75)
			pamk(sort(pe),2)->clus
			sort(log(pe))[which(clus[[1]]$partitioning==2)[1]-1]->threshold
			abline(h=threshold-shift,lty=2)
		return(list(tab,clus[[1]]$clusinfo))
	}
	else{return(tab)}
}	

