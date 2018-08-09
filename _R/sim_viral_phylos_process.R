#############################################################################################
######EXPANSIVENESS PROPERTIES ON TREES UNDER DIFFERENT TIME-DEP PROCESSES#######################
#############################################################################################


setwd('~/Desktop/_mtree/_simulated/_proc')
library(igraph)
source('sim_sample_tree.R')
source('~/Documents/GitHub/pflam/_R/_JSD.R')
setwd('~/Desktop/_mtree/_simulated/_proc/_ultrametric')

#define constant, linear, and exponential speciation functions
l1<-function(t,y){y[1]}
l3<-function(t,y){y[1]*exp(y[2]*t)}

#define constant, linear, and exponential extinction functions
m0<-function(t,y){0}

#define speciation parameters
lp1<-c(0.01)
lp6<-c(0.01,0.2)
lp7<-c(0.01,-0.1)

#define extinction parameters
mp0<-c(0)

#set tips and ages
time<-c(100,100,100)
tips<-c(50,50,50)
samples=3

#SPECIATION-EXTINCTION SCENARIOS


#stable constant birth-death
f.lambs1<-list(l1,l1,l1);lamb_pars1<-list(lp1,lp1,lp1)
#stable positive exponential time-dep (nearly) pure-birth
f.lambs2<-list(l3,l3,l3);lamb_pars2<-list(lp6,lp6,lp6)
#stable negative exponential time-dep (nearly) pure-birth
f.lambs3<-list(l3,l3,l3);lamb_pars3<-list(lp7,lp7,lp7)
	
f.lambs<-list(f.lambs1,f.lambs2,f.lambs3)
lamb_pars<-list(lamb_pars1,lamb_pars2,lamb_pars3)	

#no extinction	
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

c()->SIMs
for(j in 1:length(f.lambs)){
SIMs[[j]]<-lapply(1:100,function(q){
	print(q)
	sim_sample_tree(f.lambs[[j]],f.mus,lamb_pars[[j]],mu_pars,
		samples,time,tips,ultrametric=T)
	})				
}		


SIM<-c(SIMs[[1]], SIMs[[2]], SIMs[[3]])
names<-c(rep('CBD',100),rep('ETDp',100),rep('ETDn',100))

e<-lapply(SIM,gete)
c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(e)){
	print(n)
	#dens(log(e[[n]]+1e-9))->d[[n]]
	dens(e[[n]])->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}

	
tab<-cbind(log(pe),skew,height,names)	
colnames(tab)<-c('pe','skew','height','proc')


pdf('sim_viral_phylos_proc.pdf',12,6)
par(mfcol=c(1,3),mar=c(1,1,1,1))
lapply(seq(1,length(SIMs),1),function(p){plot(ladderize(SIMs[[p]][[p]]),show.tip.label=F)})
dev.off()

setwd('~/Desktop/_mtree/_simulated/_proc/_process(lambdastar)')
write.table(tab,file="sim_viral_phylos_table_ultrametric_proc.txt")


tabs<-read.table('sim_viral_phylos_table_ultrametric_proc.txt',header=T,row.names=1)
cols<-c('#207F71','#BB5CA3','#C2D24A')

pdf('sim_viral_phylos_proc_phylospace.pdf')
sactter.grid(tabs$pe,tabs$skew,tabs$skew*2.103423,color=cols[as.numeric(tabs$proc)],pch=20,angle=130)
	#legend("topleft",levels(tabs$proc),col=cols,pch=20,bty="n")	
dev.off()

subs<-lapply(1:3,function(j){
	subset(tabs,tabs$proc==levels(tabs$proc)[j])
})

pdf('sim_viral_phylos_proc_phylostats.pdf',6,12)
par(mfrow=c(3,1))
c()->ps
pes<-lapply(1:length(subs),function(p){c(ps,subs[[p]]$pe)->ps})
	boxplot(pes,names=levels(tabs$proc),las=2,outline=F,col=cols)
c()->sk
skews<-lapply(1:length(subs),function(p){c(sk,subs[[p]]$skew)->sk})
	boxplot(skews,names=levels(tabs$proc),las=2,outline=F,col=cols)	
c()->hs
heights <-lapply(1:length(subs),function(p){c(hs,subs[[p]]$skew*2.103423)->hs})
	boxplot(heights,names=levels(tabs$proc),las=2,outline=F,col=cols)
dev.off()	
