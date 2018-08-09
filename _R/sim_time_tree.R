require(phytools)
	require(geiger)
	require(TreeSim)

sim_sample_tree<-function(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,root,ultrametric=F,rescale=F){
	#simulate three ultrametric trees 
	if(ultrametric==T){
		trees<-lapply(1:samples,function(s){
			sim.bd.taxa(tips[s],1,f.lambs[[s]](time[s],lamb_pars[[s]]),f.mus[[s]](time[s],mu_pars[[s]]),complete=F)[[1]]
		})
	}
	#simulate three non-ultrametric trees
	else{
		tree<-lapply(1:samples,function(s){
			sim.bd.taxa(tips[s],1,f.lambs[[s]](time[s],lamb_pars[[s]]),f.mus[[s]](time[s],mu_pars[[s]]),complete=T)[[1]]
		})
		trees<-lapply(1:samples,function(t){
			drop.tip(tree[[t]],sample(abs(tips[t]-1-length(tree[[t]]$tip.label))))
		})
	}
	if(rescale==T){
	#scale trees
	trees<-lapply(1:samples,function(t){
		trees[[t]]$edge.length<-trees[[t]]$edge.length/max(nodeHeights(trees[[t]])[,2])*time[t]
		trees[[t]]
	})
	}
	#graph trees sequentially to terminal edges
	bind.tree(trees[[1]],trees[[2]],where=max(trees[[1]]$edge))->tr1
	bind.tree(tr1,trees[[3]],where=max(trees[[2]]$edge))->tr2
	return(tr2)	
}		

sim_sample_tree_root<-function(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,root,ultrametric=F,rescale=F){
	#simulate three ultrametric trees 
	if(ultrametric==T){
		trees<-lapply(1:samples,function(s){
			sim.bd.taxa(tips[s],1,f.lambs[[s]](time[s],lamb_pars[[s]]),f.mus[[s]](time[s],mu_pars[[s]]),complete=F)[[1]]->tr
			tr$root.edge<-root[s]
			tr
		})
	}
	#simulate three non-ultrametric trees
	else{
		tree<-lapply(1:samples,function(s){
			sim.bd.taxa(tips[s],1,f.lambs[[s]](time[s],lamb_pars[[s]]),f.mus[[s]](time[s],mu_pars[[s]]),complete=T)[[1]]->tr
			tr$root.edge<-root[s]
			tr
		})
		trees<-lapply(1:samples,function(t){
			drop.tip(tree[[t]],sample(abs(tips[t]-1-length(tree[[t]]$tip.label))))
		})
	}
	if(rescale==T){
	#scale trees
	trees<-lapply(1:samples,function(t){
		trees[[t]]$edge.length<-trees[[t]]$edge.length/max(nodeHeights(trees[[t]])[,2])*time[t]
		trees[[t]]
	})
	}
	#graph trees sequentially to terminal edges
	bind.tree(trees[[1]],trees[[2]],where=max(trees[[1]]$edge))->tr1
	bind.tree(tr1,trees[[3]],where=max(trees[[2]]$edge))->tr2
	return(tr2)	
}		


gete<-function(phy){
abs(eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=F),
	only.values=T)$values)
}

getn<-function(phy){
abs(eigen(
	graph.laplacian(
		graph.adjacency(
			data.matrix(dist.nodes(phy)),
			weighted=T),
		normalized=T),
	only.values=T)$values)
}

sim_sample_tree2<-function(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=F){
	require(phytools)
	require(geiger)
	require(TESS)
	if(ultrametric==F){
	trees<-lapply(1:samples,function(s){
		sim_time_tree(f.lambs[[s]],f.mus[[s]],lamb_pars[[s]],mu_pars[[s]],
			time.stop=time[s])->tr
			if(length(tr$tip.label)==0){break}
			})
		}	
	else{
		trees<-lapply(1:samples,function(s){
			tess.sim.taxa.age(1,tips,time[s],f.lambs[[s]](time[s],lamb_pars[[s]]),f.mus[[s]](time[s],mu_pars[[s]]))[[1]]
			})		
		}
	bind.tree(trees[[1]],trees[[2]],where=max(trees[[1]]$edge))->tr
	bind.tree(tr,trees[[3]],where=max(tr$edge))	
	}

sim_time_tree<-function(f.lamb,f.mu,lamb_par,mu_par,time.stop=tot_time,return.extinct=T,prune.extinct=F){
	while(1){
			nblineages<-c(1)
			times<-c(0)
			b<-f.lamb(0,lamb_par)
			d<-f.mu(0,mu_par)
			dt<-rexp(1,(b+d))
			t<-dt
			if(t>=time.stop){
				t<-time.stop
				alive<-1
				times<-c(times,t)
				nblineages<-c(nblineages,1)
				break
			}
			r<-runif(1)
			if(r>b/(b+d)){
				times<-c(times,dt)
				nblineages<-c(nblineages,0)
				alive<-rep(FALSE,1)
			}
			else{
				edge<-rbind(c(1,2),c(1,3))
				edge.length<-rep(NA,2)
				stem.depth<-rep(t,2)
				alive<-rep(TRUE,2)
				times<-c(times,dt)
				nblineages<-c(nblineages,sum(alive))
				next.node<-4
				repeat{
					if(sum(alive)==0)
					break
					b<f.lamb(t,lamb_par)
					d<-f.mu(t,mu_par)
					dt<-rexp(1,sum(alive)*(b+d))
					t<-t+dt
					if(t>=time.stop){
						t<-time.stop
						times<-c(times,t)
						nblineages<-c(nblineages,sum(alive))
						break
				}
			r<-runif(1)
			if(r<=b/(b+d)){
				random_lineage<-round(runif(1,1,sum(alive)))
				e<-matrix(edge[alive,],ncol=2)
				parent<-e[random_lineage,2]
				alive[alive][random_lineage]<-FALSE
				edge<-rbind(edge,c(parent,next.node),c(parent,next.node+1))
				next.node<-next.node+2
				alive<-c(alive,TRUE,TRUE)
				stem.depth<-c(stem.depth,t,t)
				x<-which(edge[,2]==parent)
				edge.length[x]<-t-stem.depth[x]
				edge.length<-c(edge.length,NA,NA)
				times<-c(times,t)
				nblineages<-c(nblineages,sum(alive))
			}	
			else{
				random_lineage<-round(runif(1,1,sum(alive)))
				edge.length[alive][random_lineage]<-t-stem.depth[alive][random_lineage]
				alive[alive][random_lineage]<-FALSE
				times<-c(times,t)
				nblineages<-c(nblineages,sum(alive))
				}
			}
		}	
			if(return.extinct==TRUE | sum(alive)>0){
				break
			}
		}	
			if(sum(alive)==0) {obj<-NULL}
			else if (sum(alive)==1)	{obj<-1}
			else{
				edge.length[alive]<-t-stem.depth[alive]
				n<- -1
				for(i in 1:max(edge)){
					if(any(edge[,1]==i)){
						edge[which(edge[,1]==i),1]<-n
						edge[which(edge[,2]==i),2]<-n
						n<-n-1
					}
				}
				edge[edge>0]<-1:sum(edge>0)
				tip.label<-1:sum(edge>0)
				mode(edge)<-'character'
				mode(tip.label)<-'character'
				obj<-list(edge=edge,edge.length=edge.length,tip.label=tip.label)
				class(obj)<-'phylo'
				obj<-old2new.phylo(obj)
				if(prune.extinct){
					obj<-drop.extinct(obj)
				}
			}
			return(obj)
}
			

###EXAMPLE
#library(geiger)
#l1<-function(t,y){y[1]}
#l2<-function(t,y){y[1]*exp(y[2]*t)}
#f.lambs<-list(l1,l2,l2)
#lpar1<-c(0.2)
#lpar2<-c(0.2,0.1)
#lamb_pars<-list(lpar1,lpar2,lpar2)

#m1<-function(t,y){0}
#m2<-function(t,y){0.01}
#f.mus<-list(m1,m2,m1)
#mpar1<-c()
#mu_pars<-list(mpar1,mpar1,mpar1)

#samples<-3
#time<-c(10,10,10)
#tips=10

#sim_sample_tree(f.lambs,f.mus,lamb_pars,mu_pars,samples,time,tips,ultrametric=F)->tte			
