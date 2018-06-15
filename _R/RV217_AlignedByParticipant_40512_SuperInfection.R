library(RPANDA)
library(igraph)
library(phytools)
library(geiger)
library(RColorBrewer)
source('~/Desktop/_mtree/_R/_JSD.R', chdir = TRUE)
setwd('~/Desktop/_mtree/RV217/AlignedByParticipant/40512_superinfection')

tr<-read.nexus('~/Desktop/_mtree/RV217/AlignedByParticipant/40512_superinfection/40512allseq_rooted.nex')
drop.tip.ni(tr,tips1)->tr1
drop.tip(tr,tips1)->tr2
write.nexus(tr1,file="40512allseq_rooted_infect1.nex")
write.nexus(tr2,file="40512allseq_rooted_infect2.nex")