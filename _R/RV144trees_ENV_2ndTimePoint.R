tr<-read.tree('~/Desktop/_mtree/RV144/RV144_ENV_2ndTimePoint.1626seq.fas.treefile')
inf<-read.table('~/Desktop/_mtree/RV144/RV144.collectiondates.2ndTimePoint.txt',header=T)

RV144trees<-lapply(1:dim(inf)[1],function(i){
	print(i)
	grep(inf[i,1],tr$tip.label)->tips
	drop.tip.ni(tr,tr$tip.label[tips])
	})
save(RV144trees,file="RV144trees_ENV_2ndTimePoint.RData")