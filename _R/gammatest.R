library(RPANDA)
library(phytools)
library(phylosim);PSIM_FAST <- TRUE
setwd('~/Desktop/_mtree/_simulated')
source('seqsim.R')
load('RV217bestMCCtrees_genes.RData')
setwd('~/Desktop/_mtree/_simulated/_gammatest')
genetree<-rep(c(genetrees[[1]],genetrees[[3]]),50)
#################simulate_alignments_under_different_gamma_profiles##################################simulate_alignments_under_different_gamma_profiles##################################simulate_alignments_under_different_gamma_profiles##################################simulate_alignments_under_different_gamma_profiles


gtrNormal<-lapply(1:100,function(l){
	seqsim(genetree[[l]],length=100,model='GTR',option="normal",basefreq=rep(0.25,4))
	})

gtrGamma<-lapply(1:100,function(l){
	print(l)
	seqsim(genetree[[l]],length=100,model='GTR',option="discrete",basefreq=rep(0.25,4))
	})

gtrIGamma<-lapply(1:100,function(l){
	print(l)
	seqsim(genetree[[l]],length=100,model='GTR',option="invariant",basefreq=rep(0.25,4))
	})		

gtrs<-c(gtrNormal,gtrGamma,gtrIGamma)


lapply(1:length(gtrs),function(i){
    saveAlignment(gtrs[[i]],file=paste("GTRs",i,".fas",sep=""),skip.internal=T)
  })


#################align_and_build_trees##################################align_and_build_trees##################################align_and_build_trees##################################align_and_build_trees##########

#align with mafft
for fasta_file in $(ls *.fas)
do
linsi $fasta_file > $fasta_file.out
done

#build trees
for msa_file in $(ls *.fas.out)
do
./iqtree -s $msa_file
done


#################compute_spectral_density_profile##################################compute_spectral_density_profile##################################compute_spectral_density_profile################################

#read trees
files<-mixedsort(list.files(path="~/Desktop/_mtree/_simulated/_gammatest",pattern="*.treefile"))
trees<-lapply(files,read.tree)

##compute spectral density profile
library(igraph)
source('~/Desktop/_mtree/_R/_JSD.R', chdir = TRUE)

#set SDP function
gete<-function(phy){
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

e<-lapply(trees,gete)

c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(e)){
	dens(log(e[[n]]))->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}










