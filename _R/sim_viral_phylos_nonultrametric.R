setwd('~/Desktop/_simulated/_proc')
library(igraph)
source('sim_time_tree.R')
source('~/Documents/GitHub/pflam/_R/_JSD.R')
setwd('~/Desktop/_simulated/_proc/_nonultrametric')

#define constant, linear, and exponential speciation functions
l1<-function(t,y){y[1]}
l2<-function(t,y){y[1]*t}
l3<-function(t,y){y[1]*exp(y[2]*t)}

#define constant, linear, and exponential extinction functions
m1<-function(t,y){y[1]}
m2<-function(t,y){y[1]*t}
m3<-function(t,y){y[1]*exp(y[2]*t)}

#define speciation parameters
lp1<-c(0.01)
lp2<-c(0.2)
lp3<-c(0.5)
lp4<-c(0.05,0.05)
lp5<-c(0.01,0.5)
lp6<-c(0.01,1)


#define extinction parameters
mp1<-c(0.005)
mp2<-c(0.075)
mp3<-c(0.1)
mp4<-c(0.01,0.1)
mp5<-c(0.01,0.3)
mp6<-c(0.01,-0.01)
mp7<-c(0.01,-0.05)

#set tips and ages
time<-c(10,10,10)
tips<-c(10,10,10)
samples=3

#SPECIATION-EXTINCTION SCENARIOS

##CBD
#stable constant birth-death
f.lambs<-list(l1,l1,l1)
	lamb_pars<-list(lp1,lp1,lp1)
f.mus<-list(m1,m1,m1)
	mu_pars<-list(mp1,mp1,mp1)

CBDs<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=F)
	})				
		

##ETD
#stable positive exponential time-dep (nearly) pure-birth
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp4,lp4,lp4)
f.mus<-list(m1,m1,m1)
	mu_pars<-list(mp2,mp2,mp2)

ETDs<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=F)
	})					

#increasing exponential time-dep birth-death
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp4,lp5,lp6)
f.mus<-list(m1,m1,m1)
	mu_pars<-list(mp2,mp2,mp2)

ETDi<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=F)
	})	

#decreasing exponential time-dep birth-death
f.lambs<-list(l3,l3,l3)
	lamb_pars<-list(lp6,lp5,lp4)
f.mus<-list(m3,m3,m3)
	mu_pars<-list(mp4,mp4,mp4)

ETDd<-lapply(1:200,function(q){
	print(q)
	sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=F)
	})					

simtrees<-c(CBDs,ETDs,ETDi,ETDd)

simnames<-c(rep("CBDs",200),rep("ETDs",200),rep("ETDi",200),rep("ETDd",200))#,rep("Ocle",200),rep("Oelc",200),rep("Ocbped",200),rep("Ocbned",200))

pdf('simtrees_nonultrametric.pdf',10,5)
par(mfcol=c(1,4),mar=c(1,1,1,1))
lapply(seq(1,800,200),function(p){plot(ladderize(simtrees[[p]]),show.tip.label=F,main=paste(simnames[p]))})
dev.off()

gete<-function(phy){
abs(eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=F),
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

	
tab<-cbind(log(pe),skew,height,simnames)	
colnames(tab)<-c('pe','skew','height','proc')
write.table(tab,file="sim_viral_phylos_table_nonultrametric.txt")


	
