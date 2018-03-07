get_extinction_list=function(N,tree,equal.ab=TRUE){
  
  tip.ab=tree$tip.ab
  names(tip.ab) = tree$tip.label
  tips.order = NULL
  if (equal.ab == TRUE) {
    tips.order = sample(names(tip.ab), N)
    
  } else {  
    ab_null = c(which(is.na(tip.ab)), which(tip.ab==0))
    if (length(ab_null)>0) {
      ab_null_order = sample(ab_null, length(ab_null))
      tip_not_null = tip.ab[-ab_null]
      tips.order = c(names(ab_null_order), names(tip_not_null)[order(tip_not_null,names(tip_not_null))]) 
    } else {
      tips.order = names(tip.ab)[order(tip.ab,names(tip.ab))] 
    }
  } 
  return(tips.order)
}

get_PD_sample=function(epsilon,beta,alpha,N,sampl.frac,ntree,equal.ab,eta,lengths="yule",b=1,d=0){
  lambda=lambda_N(epsilon,beta,N,upper=10000)
  PD=pbsapply(1:ntree,function(i){
    if(lengths=="yule"){
      tree=simulate_yule(epsilon,alpha,beta,N,b = b,d=d,equal.ab = equal.ab,eta = eta,lambda=lambda)
    }else if(lengths=="kingman"){
      tree=simulate_kingman(epsilon,alpha,beta,N,equal.ab = equal.ab,eta = eta,lambda=lambda)
    }
    M=max(depth(tree))
    tip=get_extinction_list(N,tree,equal.ab)
    PD.ech=sapply(1:length(sampl.frac),function(j){
      if (round(sampl.frac[j]*N)==N){sum(depth(tree))+M
      } else if (round(sampl.frac[j]*N) == 0) {
        0
      } else if (round(sampl.frac[j]*N) == 1) {
        M
      } else {
        nd=depth(drop.tip(tree,tip[1:(N-round(sampl.frac[j]*N))]))
        sum(nd)+M}})
    return(PD.ech/(sum(depth(tree))+M))
  })
  rownames(PD) = sampl.frac
  colnames(PD) = seq(1,ntree)
  return(PD)
}
