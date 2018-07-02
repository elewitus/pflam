#mat=n x m matrix with column names
#no NAs (convert to 0?)
#unsupervised clustering for range of 1:Y clusters 
#using expectation maximization on different gaussian mixture models
#outputs plots of expression through time for each cluster and cohesion coefficient for each gene
#outputs BIC scores for clusters 1:Y (under best fit model)
#outputs dataframes of clusters

brassica<-function(mat,Y){
	#get CI function
	CIs<-function(mean,sd,size){
		rbind(mean+qnorm(0.975)*sd/sqrt(size),
			mean-qnorm(0.975)*sd/sqrt(size))
	}
	#call required packages
	require(gtools)
	require(cluster)
	require(mclust)
	#get BIC scores for 1:Y clusters under spherical,ellipsoidal,diagonal models
	Mclust(mat,G=1:Y)->j
	print(paste("best model fit",j[[3]],"with",j[[6]],"clusters"))
	#assign genes to clusters and calculate silhouette widths for each gene
	pam(mat,j[[6]])->getClust
	getClust$silinfo$widths[mixedsort(rownames(getClust$silinfo$widths)),]->rtab
	cbind(mat,rtab[,c(1,3)])->newmat
	#subset mat by cluster
	newmats<-lapply(1:j[[6]],function(l){
		as.data.frame(subset(newmat,newmat[,dim(mat)[2]+1]==l))->newmats
		colnames(newmats)<-c(colnames(mat),"cluster","silinfo")
		newmats
		})
	#tally median and CIs for each cluster over time
	c()->med
	for(m in 1:length(newmats)){	
		med[[m]]<-sapply(1:dim(mat)[2],function(n){
			median(newmats[[m]][,n])->median
			sd<-sqrt(sum((newmats[[m]][,n]-mean(newmats[[m]][,n]))^2/(length(newmats[[m]][,n])-1)))
			rbind(median,CIs(median,sd,dim(newmats[[m]])[2]))
		})
	}
	#plot CIgons and silhouette widths for each cluster
	par(mfrow=c(round(j[[6]]),2))
		lapply(1:j[[6]],function(o){
		plot(1:dim(mat)[2],med[[o]][1,],type='l',lwd=2,
			ylim=c(range(med[[o]])),ann=F,axes=F)
			polygon(c(1:dim(mat)[2],rev(1:dim(mat)[2])),
					c(med[[o]][3,],rev(med[[o]][2,])),
					col="light grey")
			axis(1,at=1:dim(mat)[2],labels=colnames(mat))	
			axis(2,las=2)
			mtext('reads',2,2)	
		barplot(sort(newmats[[o]]$silinfo),
			horiz=T,col=NA,main=paste(dim(newmats[[o]])[1],"genes in cluster; average width=",round(mean(newmats[[o]]$silinfo),2)))		
			})
			mtext("silhouette width",1,2.5)
	res<-list(as.table(j[[7]])[,which(colnames(as.table(j[[7]]))==j[[3]])],newmats)
	#res[[1]]: BIC scores for 1:Y clusters under best fit model
	#res[[2]][[n]]: dataframes of gene clusters 1:Y
	return(res)		
}



