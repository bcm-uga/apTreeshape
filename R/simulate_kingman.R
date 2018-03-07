build_tree=function(H,tip.lab=rep(NA,length(H)+1)){
  if (length(H)==0) {
    tree=list(edge=matrix(c(2,1),1,2),edge.length=0,Nnode=1,tip.label=tip.lab)
    class(tree)="phylo"
    return(tree)
  }else{
    H.crown=max(H)
    split=which.max(H)
    t1=tip.lab[1:split]
    t2=tip.lab[(split+1):length(tip.lab)]
    if (split>1) {H1=H[1:(split-1)]
    }else{H1=c()}
    if (split<length(H)){H2=H[(split+1):length(H)]
    }else {H2=c()}
    
    a1=build_tree(H1,t1)
    a2=build_tree(H2,t2)
    root1=H.crown-max(c(H1,0))
    root2=H.crown-max(c(H2,0))
    return(collapse.singles(bind.trees(a1,a2,root1,root2)))
  }
}

depth=function(tree){
  d=node.depth.edgelength(tree)
  d[1]-d[(tree$Nnode+2):(2*tree$Nnode+1)]
}

nodes_depths_ordonnes=function(tree){
  aux=function(phylo,noeud=phylo$Nnode+2){
    if(phylo$Nnode<=1){return(max(phylo$edge.length))
    }else{
      fils=phylo$edge[phylo$edge[,1]==noeud,][,2]
      if (fils[1]<(phylo$Nnode+2)){
        ad=extract.clade(phylo,fils[2])
        return(cbind(max(phylo$edge.length[phylo$edge[,1]==noeud]),aux(ad)))
      }else if (fils[2]<(phylo$Nnode+2)){
        ag=extract.clade(phylo,fils[1])
        return(cbind(aux(ag),max(phylo$edge.length[phylo$edge[,1]==noeud])))
      }else{
        ag=extract.clade(phylo,fils[1])
        ad=extract.clade(phylo,fils[2])
        hg=aux(ag)
        hd=aux(ad)
        h=max(hg)+phylo$edge.length[phylo$edge[,2]==fils[1]]
        return(cbind(hg,h,hd))
      }
    }
  }
  H=aux(tree)
  return(H)
}

simulate_kingman=function(epsilon,alpha,beta,N,equal.ab=TRUE,eta=1,lambda=NULL){
  tree=simulate_tree(epsilon,alpha,beta,N,equal.ab,eta,lambda=lambda)
  ab=tree$tip.ab
  order=rank(nodes_depths_ordonnes(tree))
  depths=cumsum(rexp(N-1,sapply((N:2),function(i){i*(i-1)/2})))
  tree=build_tree(depths[order],1:N)
  tree$tip.ab=ab
  return(tree)
}
