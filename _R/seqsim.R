codsim<-function(tree,length,titv,dnds_min,dnds_max,plot=T){
	p<-GY94()
	p$kappa<-titv
	codon.freqs <- abs(rnorm(61, mean = 10, sd = 3))
	codon.freqs <- codon.freqs/sum(codon.freqs)
	p$equDist <- codon.freqs
	s<-CodonSequence(length=length,processes=list(list(p)))
	sampleStates(s)
	omegaVarM3(s,p,omegas=c(runif(4,dnds_min,dnds_max)),probs=rep(0.25,4))
	if(plot==T){
		omegaHist(s,p)}
	sim <- PhyloSim(root.seq = s,phylo=tree)
	Simulate(sim)	
}

seqsim<-function(tree,length,basefreq,alpha,beta,model=c(''),option=c(''),superinfect=F,plot=F){
	options(warn=-1)
	tree$tip.label<-1:length(tree$tip.label)

	##JC69
	##Equal substitution rates and equal base frequencies (Jukes and Cantor, 1969)
	if(model=='JC69'){
		mod<-JC69()
	}
	##HKY 
	##Unequal transition/transversion rates and unequal base freq. (Hasegawa, Kishino and Yano, 1985)
	if(model=='HKY'){
		mod<-HKY(rate.params=list("Alpha"=alpha,"Beta"=beta), base.freqs=basefreq)
			}
	##F81
	##Equal rates but unequal base freq. (Felsenstein, 1981)
	if(model=='F81'){
		mod<-list(F81(base.freqs=basefreq))
			}
	##K80
	##Unequal transition/transversion rates and equal base freq. (Kimura, 1980)
	if(model=='K80'){
		mod<-K80(rate.params=list("Alpha"=alpha,"Beta"=beta), base.freqs=basefreq)
			}
	##GTR
	##General time reversible model with unequal rates and unequal base freq. (Tavare, 1986)
	if(model=='GTR'){
		mod<-GTR(rate.params=list('a'=1,'b'=2,'c'=3,'d'=1,'e'=2,'f'=3),base.freqs=basefreq)
			}
	if(option=='discrete'){
		sim<-PhyloSim(sampleStates(plusGamma(
					NucleotideSequence(len=length,proc=list(list(mod))),
					mod,0.5)),phylo=tree)
			}
	if(option=='invariant'){
		sim<-PhyloSim(sampleStates(plusInvGamma(
					NucleotideSequence(len=length,proc=list(list(mod))),
					mod,pinv=0.8,0.5)),phylo=tree)
			}
	if(option=='normal'){
		sim<-PhyloSim(sampleStates(NucleotideSequence(len=length,proc=list(list(mod)))),
			phylo=tree)
			}
	if(superinfect==T){
			super<-function(seq){
				if(!isAttached(seq$sites[[1]],mod)){
					return(seq);
				}
				cat("Resampling rate multipliers");
				plusGamma(seq,mod,1)
				return(seq)
				}
			attachHookToNode(sim,node=length(tree$tip.label)+1,fun=super)
		alignment<-Simulate(sim)
			}	
	if(superinfect==F){
		alignment<-Simulate(sim)
			}
	if(plot==T){
		pdf('seqsim.pdf')
		par(mfrow=c(1,2))
		plot(mod)
		plot(sim,num.pages=1)
		dev.off()
	}		
	options(warn=0)		
	return(alignment)	
} 	
