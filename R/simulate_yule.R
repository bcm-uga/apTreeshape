yule_lengths=function(N,b,d,tmax){
  r=b-d
  
  f=function(u){
    if (r==0){
      u/(b*(1-u))
    }else{
      log((r*u/(b*(1-u)))+1)/r
    }
  }
  
  X=c(tmax)
  n=1
  while (n<N){
    t=f(runif(1,0,1))
    if (t<tmax){X=c(X,t);n=n+1}
  }
  
  return(X[2:N])
}

yule=function(N,b,d,tmax){
  build_tree(yule_lengths(N,b,d,tmax))
}

simulate_yule=function(epsilon,alpha,beta,N,b,d,tmax=Inf,equal.ab=TRUE,eta=1,lambda=NULL){
  tree=simulate_tree(epsilon,alpha,beta,N,equal.ab,eta,lambda=lambda)
  ab=tree$tip.ab
  order=rank(nodes_depths_ordonnes(tree))
  depths=sort(yule_lengths(N,b,d,tmax))
  tree=build_tree(depths[order],1:N)
  tree$tip.ab=ab
  return(tree)
}
