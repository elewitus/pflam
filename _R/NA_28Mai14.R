tr<-read.tree('~/Desktop/_mtree/RV254/RV254.NA.28mai14.FASTA.treefile')
inf<-read.table('~/Desktop/_mtree/RV254/RV254_info.txt',header=T)

RV254trees<-lapply(1:dim(inf)[1],function(i){
	print(i)
	grep(inf[i,1],tr$tip.label)->tips
	drop.tip.ni(tr,tr$tip.label[tips])
	})
save(RV254trees,file="RV254trees_NA_28mai14.RData")
