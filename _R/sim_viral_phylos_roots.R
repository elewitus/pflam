#############################################################################################
######BETWEEN/WITHIN VARIANCE ON CBD TREES WITH DIFFERENT ROOT LENGTHS#######################
#############################################################################################


setwd('~/Desktop/_mtree/_simulated/_proc')
library(igraph)
source('sim_time_tree.R')
source('~/Documents/GitHub/pflam/_R/_JSD.R')
setwd('~/Desktop/_mtree/_simulated/_proc/_ultrametric')

#define constant, linear, and exponential speciation functions
l1<-function(t,y){y[1]}

#define constant, linear, and exponential extinction functions
m0<-function(t,y){0}

#define speciation parameters
lp1<-c(0.01)

#define extinction parameters
mp0<-c(0)

#set tips and ages
time<-rep(100,3)
tips<-rep(50,3)
samples=3

#SPECIATION-EXTINCTION SCENARIOS

#average between variance
root1<-rep(0,3)
root2<-rep(50,3)
root3<-rep(200,3)
root4<-rep(500,3)
roots<-list(root1,root2,root3,root4)

##CBD
#stable constant birth-death
f.lambs<-list(l1,l1,l1)
	lamb_pars<-list(lp1,lp1,lp1)
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

c()->CBDs
for(r in 1:length(roots)){
CBDs[[r]]<-lapply(1:100,function(q){
	print(q)
	sim_sample_tree_root(f.lambs,f.mus,lamb_pars,mu_pars,
		samples,time,tips,roots[[r]],ultrametric=T)
	})				
}	

CBD<-c(CBDs[[1]],CBDs[[2]],CBDs[[3]],CBDs[[4]])
names<-c(rep('1',100),rep('2',100),rep('3',100),rep('4',100))

e<-lapply(CBD,gete)
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


pdf('sim_viral_phylos_CBD_roots.pdf',12,6)
par(mfcol=c(1,4),mar=c(1,1,1,1))
lapply(seq(1,length(CBDs),1),function(p){plot(ladderize(CBDs[[p]][[p]]),show.tip.label=F)})
dev.off()

setwd('~/Desktop/_mtree/_simulated/_proc/_rootVariance(eta)')
write.table(tab,file="sim_viral_phylos_table_CBD_ultrametric_roots.txt")


tabs<-read.table('sim_viral_phylos_table_CBD_ultrametric_roots.txt',header=T,row.names=1)
cols<-seq(15,60,15)

pdf('sim_viral_phylos_root_phylospace.pdf')
sactter.grid(tabs$pe,tabs$skew,tabs$height,color=colors(1)[as.numeric(tabs$proc)*15],pch=20,angle=30)
	#legend("topleft",c('1','2','3','4'),col=colors(1)[seq(15,60,15)],pch=20,bty="n")	
dev.off()

subs<-lapply(1:4,function(j){
	subset(tabs,tabs$proc==j)
})

pdf('sim_viral_phylos_CBD_roots_phylostats.pdf',6,12)
par(mfrow=c(3,1))
c()->ps
pes<-lapply(1:length(subs),function(p){c(ps,subs[[p]]$pe)->ps})
	boxplot(pes,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])
c()->sk
skews<-lapply(1:length(subs),function(p){c(sk,subs[[p]]$skew)->sk})
	boxplot(skews,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])	
c()->hs
heights <-lapply(1:length(subs),function(p){c(hs,exp(subs[[p]]$height))->hs})
	boxplot(heights,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])
dev.off()	
