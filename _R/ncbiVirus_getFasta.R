library(ape)
library(seqinr)
library(rentrez)
setwd('~/Desktop/_newClades/_virus')
vf<-read.table('ncbi_virus_families.txt')
vfr<-read.table('vfr',row.names=1)

lapply(1:dim(vf)[1],function(o){
  print(o)
	entrez_search(db="nuccore",term=vfr[o,],retmax=500)->search
  if(length(search$ids)>20){
	entrez_fetch(db="nuccore", id=search$ids, rettype="fasta")->seqs
	write(seqs,file=paste(vf[o,],".fasta"),sep='\n')
  }
  else{}
})

#rename
#for f in *\ *; do mv "$f" "${f// /_}"; done
