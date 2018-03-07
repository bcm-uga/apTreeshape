get_tree_beta=function(epsilon,beta,alpha,N,sampl.frac,ntree,equal.ab=TRUE,eta=1){
  beta=pbsapply(1:ntree,function(i){
    tree=simulate_tree(epsilon, alpha, beta, N, equal.ab=equal.ab, eta=eta)
    tree$tip.label=1:length(tree$tip.label)
    tip=get_extinction_list(N,tree,equal.ab)
    return(sapply(1:length(sampl.frac),function(j){
      if (round(sampl.frac[j]*N)==N){maxlik.betasplit(as.treeshape(tree), up=10, remove.outgroup=FALSE, confidence.interval="none")$max_lik
      } else if (round(sampl.frac[j]*N)==0 || round(sampl.frac[j]*N)==1) {
        NA
      } else {
        maxlik.betasplit(as.treeshape(drop.tip(tree,tip[1:(N-round(sampl.frac[j]*N))])), up=10, remove.outgroup=FALSE, confidence.interval="none")$max_lik       
      }
    }))
  })
  rownames(beta) = sampl.frac
  colnames(beta) = seq(1,ntree)
  return(beta)
}
