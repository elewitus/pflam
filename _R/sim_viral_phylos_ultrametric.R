setwd('~/Desktop/_simulated/_proc')
library(igraph)
source('sim_time_tree.R')
source('~/Documents/GitHub/pflam/_R/_JSD.R')

#define constant, linear, and exponential speciation functions
l1<-function(t,y){y[1]}
l3<-function(t,y){y[1]*exp(y[2]*t)}

#define constant, linear, and exponential extinction functions
m0<-function(t,y){0}
m1<-function(t,y){y[1]}
m3<-function(t,y){y[1]*exp(y[2]*t)}

#define speciation parameters
lp1<-c(0.01)
lp2<-c(0.2)
lp3<-c(0.5)
lp4<-c(0.05,0.05)
lp5<-c(0.01,0.5)
lp6<-c(0.01,1)
lp7<-c(0.01,-0.2)


#define extinction parameters
mp0<-c(0)
mp1<-c(0.005)
mp2<-c(0.075)
mp3<-c(0.1)
mp4<-c(0.01,0.1)
mp5<-c(0.01,0.3)
mp6<-c(0.01,-0.01)
mp7<-c(0.01,-0.05)

#set tips and ages
time<-c(100,100,100)
tips<-c(50,50,50)
samples=3

#SPECIATION-EXTINCTION SCENARIOS

##CBD
#stable constant birth-death
f.lambs<-list(l1,l1,l1)
	lamb_pars<-list(lp1,lp1,lp1)
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

CBDs<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=T)
	})				
		

##ETD
#stable positive exponential time-dep (nearly) pure-birth
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp4,lp4,lp4)
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

ETDp<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=T)
	})					

#stable negative exponential time-dep (nearly) pure-birth
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp7,lp7,lp7)
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

ETDn<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=T)
	})		
			

#increasing exponential time-dep (nearly) pure-birth
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp4,lp5,lp6)
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

ETDi<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=T)
	})	

#increasing exponential time-dep (nearly) pure-birth
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp6,lp5,lp4)
f.mus<-list(m0, m0, m0)
	mu_pars<-list(mp0, mp0, mp0)

ETDd<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=T)
	})

simtrees<-c(CBDs,ETDp,ETDn,ETDi,ETDd)

simnames<-c(rep("CBDs",200),rep("ETDp",200),rep("ETDn",200),rep("ETDi",200),rep("ETDd",200))

#pdf('simtrees_ultrametric_large.pdf',10,5)
par(mfcol=c(1,5),mar=c(1,1,1,1))
lapply(seq(1,length(simtrees),200),function(p){plot(ladderize(simtrees[[p]]),show.tip.label=F,main=paste(simnames[p]))})
#dev.off()

gete<-function(phy){
abs(eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=T),
	only.values=T)$values)
}

e<-lapply(simtrees,gete)

c()->d;c()->dsc;c()->pe;c()->skew;c()->height
for(n in 1:length(e)){
	dens(log(e[[n]]+1e-3))->d[[n]]
	#dens(e[[n]])->d[[n]]
	d[[n]]$y/integr(d[[n]]$x,d[[n]]$y)->dsc[[n]]
	max(e[[n]])->pe[[n]]
	skewness(dsc[[n]])->skew[[n]]
	max(dsc[[n]])->height[[n]]
	}

	
tab<-cbind(pe,skew,height,simnames)	
colnames(tab)<-c('pe','skew','height','proc')
setwd('~/Desktop/_simulated/_proc/_ultrametric')
write.table(tab,file="sim_viral_phylos_table_ultrametric_large_normal.txt")

cols<-c(1,5,4,3,2)*15
tabs<-read.table('sim_viral_phylos_table_ultrametric_large.txt',header=T,row.names=1)
sactter.grid(tabs$pe,tabs$skew,tabs$height,color=colors(1)[as.numeric(tabs$proc)*15],pch=20,angle=230)
	legend("topleft",levels(tabs$proc),col=colors(1)[cols],pch=20,bty="n")	

subs<-lapply(1:length(levels(tabs$proc)),function(j){
	subset(tabs,tabs$proc==levels(tabs$proc)[j])
})

par(mfrow=c(3,1))
c()->ps
pes<-lapply(1:length(subs),function(p){c(ps,subs[[p]]$pe)->ps})
	boxplot(pes,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])
c()->sk
skews<-lapply(1:length(subs),function(p){c(sk,subs[[p]]$skew)->sk})
	boxplot(skews,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])	
c()->hs
heights <-lapply(1:length(subs),function(p){c(hs,subs[[p]]$height)->hs})
	boxplot(heights,names=levels(tabs$proc),las=2,outline=F,col=colors(1)[cols])

par(mfrow=c(3,1))
plot(tabs$pe,tabs$height,col=colors(1)[as.numeric(tabs$proc)*15],pch=20)
plot(tabs$pe,tabs$skew,col=colors(1)[as.numeric(tabs$proc)*15],pch=20)
plot(tabs$skew,tabs$height,col=colors(1)[as.numeric(tabs$proc)*15],pch=20)