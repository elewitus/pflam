library(RPANDA)
library(phytools);PSIM_FAST <- TRUE
library(phylosim)
setwd('~/Desktop/_mtree/_simulated')
load('~/Desktop/_mtree/RV217/AlignedByParticipant/_genes/_geneBeast/RV217bestMCCtrees_genes.RData')

#gene lengths
env=2646
gag=1503
nef=633
pol=3013
rev=351
tat=306
vif=585
vpr=291
vpu=246
length<-c(rep(env,43),rep(gag,43),rep(nef,43),rep(pol,43),rep(rev,43),rep(tat,43),rep(vif,43),rep(vpr,43),rep(vpu,43))

#variants
options<-c('normal','discrete','invariant')
basefreq<-list(c(0.4,0.3,0.2,0.1),c(0.1,0.2,0.3,0.4))
alpha<-c(10,1)
beta<-c(1,10)


##JC69
c()->jc69
for(o in 1:length(options)){
	print(o)
jc69[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],model='JC69',option=options[o])	
	})
}	

lapply(1:length(jc69),function(i){
	lapply(1:length(jc69[[i]]),function(j){
		saveAlignment(jc69[[i]][[j]],file=paste("JC69",options[i],j,"fas",sep=""),skip.internal=T)
	})
})

##HKY
#basefreq1_highAlpha_lowBeta
c()->hky_basefreq1_highAlpha_lowBeta
for(o in 1:length(options)){
	print(o)
hky_basefreq1_highAlpha_lowBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='HKY',alpha=alpha[1],beta=beta[1],basefreq=basefreq[[1]],option=options[o])	
	})
#basefreq2_highAlpha_lowBeta
c()->hky_basefreq2_highAlpha_lowBeta
for(o in 1:length(options)){
	print(o)
hky_basefreq2_highAlpha_lowBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='HKY',alpha=alpha[1],beta=beta[1],basefreq=basefreq[[2]],option=options[o])	
	})
#basefreq1_lowAlpha_highBeta
c()->hky_basefreq1_lowAlpha_highBeta
for(o in 1:length(options)){
	print(o)
hky_basefreq1_lowAlpha_highBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='HKY',alpha=alpha[2],beta=beta[2],basefreq=basefreq[[1]],option=options[o])	
	})
#basefreq2_lowAlpha_highBeta
c()->hky_basefreq2_lowAlpha_highBeta
for(o in 1:length(options)){
	print(o)
hky_basefreq2_lowAlpha_highBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='HKY',alpha=alpha[2],beta=beta[2],basefreq=basefreq[[2]],option=options[o])	
	})


##F81
#basefreq1
c()->f81_basefreq1
for(o in 1:length(options)){
	print(o)
f81_basefreq1[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],model='F81',basefreq=basefreq[[1]],option=options[o])	
	})
#basefreq2
c()->f81_basefreq2
for(o in 1:length(options)){
	print(o)
f81_basefreq2[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],model='F81',basefreq=basefreq[[2]],option=options[o])	
	})	
	

##K80
#basefreq1_highAlpha_lowBeta
c()->k80_basefreq1_highAlpha_lowBeta
for(o in 1:length(options)){
	print(o)
k80_basefreq1_highAlpha_lowBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='k80',alpha=alpha[1],beta=beta[1],basefreq=basefreq[[1]],option=options[o])	
	})
#basefreq2_highAlpha_lowBeta
c()->k80_basefreq2_highAlpha_lowBeta
for(o in 1:length(options)){
	print(o)
k80_basefreq2_highAlpha_lowBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='k80',alpha=alpha[1],beta=beta[1],basefreq=basefreq[[2]],option=options[o])	
	})
#basefreq1_lowAlpha_highBeta
c()->k80_basefreq1_lowAlpha_highBeta
for(o in 1:length(options)){
	print(o)
k80_basefreq1_lowAlpha_highBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='k80',alpha=alpha[2],beta=beta[2],basefreq=basefreq[[1]],option=options[o])	
	})
#basefreq2_lowAlpha_highBeta
c()->k80_basefreq2_lowAlpha_highBeta
for(o in 1:length(options)){
	print(o)
k80_basefreq2_lowAlpha_highBeta[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],
		model='k80',alpha=alpha[2],beta=beta[2],basefreq=basefreq[[2]],option=options[o])	
	})	
	
##GTR
#basefreq1
c()->gtr_basefreq1
for(o in 1:length(options)){
	print(o)
gtr_basefreq1[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],model='gtr',basefreq=basefreq[[1]],option=options[o])	
	})
#basefreq2
c()->gtr_basefreq2
for(o in 1:length(options)){
	print(o)
gtr_basefreq2[[o]]<-lapply(1:length(genetrees),function(l){
		seqsim(genetrees[[l]],length=length[l],model='gtr',basefreq=basefreq[[2]],option=options[o])	
	})	

##
	count<-1
	trss<-list()
	for(i in 1:length(trs)){
	if(is.phylo(trs[[i]]) == TRUE){
		trs[[i]]->trss[[count]]
		count <- count + 1
	}
}

	count<-1
	list()->tr
	for(i in 1:length(trss)){
		if(length(trss[[i]]$tip.label) > 10){
			trss[[i]]->tr[[count]]
			count<-count+1
		}
	}
