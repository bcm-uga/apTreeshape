transform_eta=function(tree,unsampled,eta,list_depth=NULL,change_depth=0){
  split=split(tree)
  ntip=tree$Nnode+1
  depth=c(rep(0,ntip),depth(tree))
  A1=-1
  A2=-1
  p_unsampled=0
  compute_depth=is.null(list_depth)
  if(compute_depth){list_depth=list(c())}
  fils=tree$edge[tree$edge[,1]==(ntip+1),2]
  if(fils[1]<=ntip) A1=tree$tip.ab[fils[1]]
  if(fils[2]<=ntip) {
    if (A1==-1) {A1=tree$tip.ab[fils[2]]}else{A2=tree$tip.ab[fils[2]]}}
  M=c(ntip+1,-1,depth[ntip+1],ntip,A1,A2)
  k=2*ntip
  id.unsampled=1
  for (i in (ntip+2):(2*ntip-1)){
    A1=-1
    A2=-1
    father=tree$edge[tree$edge[,2]==i,1]
    node_father=father
    pf=depth[i]
    pp=depth[father]
    fils=tree$edge[tree$edge[,1]==i,2]
    if(fils[1]<=ntip) A1=tree$tip.ab[fils[1]]
    if(fils[2]<=ntip) {
      if (A1==-1) {A1=tree$tip.ab[fils[2]]}else{A2=tree$tip.ab[fils[2]]}}
    if (unsampled[i]>0){
      if(compute_depth | (change_depth==i)){
        depths=runif(unsampled[i],pf,pp)
        depths=depths[order(depths,decreasing = TRUE)]
        list_depth[[i]]=depths
      }else{
        depths=list_depth[[i]]
      }
      for (j in 1:unsampled[i]){
        M=rbind(M,c(k,father,depths[j],split[1,i],-1,-1))
        father=k
        k=k+1
        id.unsampled=id.unsampled+1
      }
      p_unsampled=p_unsampled+sum(log(1:unsampled[i]))-log(pp-pf)*unsampled[i]
    }
    M=rbind(M,c(i,father,pf,split[1,i],A1,A2))
  }
  
  colnames(M)=c("node","father","depth","ntip","A1","A2")
  for(j in 1:M[1,"depth"]){
    N=sum(M[,"depth"]>(j-1) & M[,"depth"]<(j))
    if(N>0){
      p_unsampled=p_unsampled-sum(log(1:N))
    }
  }
  return(list(M=M,p_unsampled=p_unsampled,list_depth=list_depth))
}

enhance_eta=function(tree,alpha,beta,eta,epsilon,lambdaN,aN,uns){
  Rmin=log(1e-10)-1
  maxChange=15
  
  N=tree$Nnode
  unsampled=c(-1,rep(-1,2*N))
  p_nUnsampled=0
  
  for (i in (N+2):(2*N+1)){
    n=extract.clade(tree,i)$Nnode+1
    
    if(n==1){unsampled[i]=-1
    }else {
      lambda=lambdaN[n]
      an=aN[n]
      sigma=rexp(1,an)
      unsampled[i]=rpois(1,lambda*sigma)
      p_nUnsampled=p_nUnsampled+log(an)-(unsampled[i]+1)*log(lambda+an)+(unsampled[i])*log(lambda)
      if(is.na(unsampled[i])) print(paste("a",alpha,"b",beta,"t",unsampled[i],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  
  tr=transform_eta(tree,unsampled,eta)
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  
  int=c(-0,rep(-Inf,nrow(tr)-1))
  Rs=c(0,rep(-Inf,nrow(tr)-1))
  As=rep(-1,nrow(tr))
  
  for(i in 1:N){
    do=(tr[,2]>=(2*N+2) & Rs==-Inf & tr[,4]==(i+1))
    if(sum(do)>0){Rs[do]=uns[[i]][sample(10000,sum(do),replace=TRUE)]}
  }
  
  nTr= nrow(tr)
  nodes=rep(TRUE,nTr)
  depth=tr[,"depth"]
  while(sum(nodes)>0){
    i=which.min(depth)
    node=tr[i,"node"]
    
    if(!((tr[i,"A1"]==-1) | (tr[i,"A2"]==-1 & tr[i,"node"]<(2*N+2)))){
      depth[i]=Inf
      if(tr[i,"A2"]>-1){
        As[i]=tr[i,"A1"]+tr[i,"A2"]
      }else{
        k=which(tr[,"father"]==node)
        As[i]=tr[i,"A1"]*(1+exp(eta*(log(1-exp(Rs[k]))-(Rs[k]))))
      }
      father=tr[i,"father"]
      nodes[i]=FALSE
      if(father>-1) {
        j=which(tr[,"node"]==father)
        if(tr[j,"A1"] ==-1){
          tr[j,"A1"] = As[i]
        }else{
          tr[j,"A2"] = As[i]
        }
      }
    }
  }
  
  if(any(is.infinite(As))){return(list(Infinite=TRUE))}
  
  
  for(i in 2:nrow(tr)){
    father=tr[i,2]
    j=which(tr[,1]==father)
    
    if(Rs[i]==-Inf){
      A=max(1e-300,(As[j]-As[i])/As[i])
      R=-log(exp((1/eta)*log(A))+1)
      R=max(R,Rmin)
      Rs[i]=R
    }
    int[i]=int[j]+Rs[i]
  }
  
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,As=As,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=FALSE))
}

change_int_eta=function(enhance,ind,tree,alpha,beta,eta,epsilon,lambdaN,aN,change=FALSE,uns){
  Rmin=log(1e-10)-1
  maxChange=15
  N=tree$Nnode
  rm=c()
  rmR=c()
  unsampled=enhance$unsampled
  int=enhance$int
  p_nUnsampled=enhance$p_nUnsampled
  list_depth=enhance$list_depth
  continue=TRUE
  if(is.null(ind)){
    rm=c()
  }else{
    n=extract.clade(tree,ind)$Nnode+1
    tr=enhance$transform
    node=which(tr[,1]==ind)
    if(unsampled[ind]>0){
      rm=(node-1):(node-unsampled[ind])
    }
    nodeR=node-length(rm)
    rmR=rm
    if(n==1){unsampled[ind]=-1
    }else if(ind==(N+2)){unsampled[ind]=-1}else{
      lambda=lambdaN[n]
      an=aN[n]
      p_nUnsampled=p_nUnsampled+(unsampled[ind]+1)*log(lambda+an)-(unsampled[ind])*log(lambda)
      u=runif(1)
      if(u<(0.4)){
        unsampled[ind]=unsampled[ind]-1-floor(rexp(1,0.2))
        ti=unsampled[ind]
      }else if(u>0.6){
        unsampled[ind]=unsampled[ind]+1+floor(rexp(1,0.2))
        ti=unsampled[ind]
      }else{
        if(unsampled[ind]>1){
          lrm=sample(min(unsampled[ind],maxChange),1)
          rmR=sample(rm,lrm)
          nodeR=node-lrm
          ti=lrm
        }else{
          ti=unsampled[ind]
        }
      }
      node=node-length(rm)
      if(unsampled[ind]<0){
        p_nUnsampled=-Inf
        continue=FALSE
      }else{
        p_nUnsampled=p_nUnsampled-(unsampled[ind]+1)*log(lambda+an)+(unsampled[ind])*log(lambda)
        continue=TRUE
      }
      
      if(is.na(unsampled[ind])) print(paste("a",alpha,"b",beta,"t",unsampled[ind],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  if(! continue){
    enhance$p_nUnsampled=-Inf
    return(enhance)
  }
  
  tr=transform_eta(tree,unsampled,eta,list_depth = list_depth,change_depth = max(ind,0))
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  
  if(length(rm)>0){
    Rs=enhance$Rs[-rmR]
    depth=enhance$transform[-rm,"depth"]
    
  }else{
    Rs=enhance$Rs
    depth=enhance$transform[,"depth"]
    
  }
  if(! is.null(ind)){
    if(unsampled[ind]>-1){
      add=rep(-Inf,ti+1)
      if(nodeR<length(Rs)){
        Rs=c(Rs[1:(nodeR-1)],add,Rs[(nodeR+1):length(Rs)])
      }else{
        Rs=c(Rs[1:(nodeR-1)],add)}
      if(node<length(depth)){
        depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"],depth[(node+1):length(depth)])
      }else{
        depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"])
      }
    }
  }
  
  tr[,"depth"]=depth
  
  Rs[tr[,"father"]<(2*N+2)]=-Inf
  if(change){Rs[]=-1}
  Rs[1]=0
  int=c(0,rep(-Inf,nrow(tr)-1))
  As=rep(-1,nrow(tr))
  
  for(i in 1:N){
    do=(tr[,2]>=(2*N+2) & Rs==-Inf & tr[,4]==(i+1))
    samp=sample(10000,sum(do),replace=TRUE)
    if(sum(do)>0){Rs[do]=uns[[i]][samp]}
  }
  
  nTr= nrow(tr)
  nodes=rep(TRUE,nTr)
  while(sum(nodes)>0){
    i=which.min(depth)
    node=tr[i,"node"]
    if((tr[i,"A1"]==-1) | (tr[i,"A2"]==-1 & tr[i,"node"]<(2*N+2))){
      print("error in change int")
    }else{
      depth[i]=Inf
      if(tr[i,"A2"]>-1){
        As[i]=tr[i,"A1"]+tr[i,"A2"]
      }else{
        k=which(tr[,"father"]==node)
        As[i]=tr[i,"A1"]*(1+exp(eta*(log(1-exp(Rs[k]))-(Rs[k]))))
      }
      father=tr[i,"father"]
      nodes[i]=FALSE
      if(father>-1) {
        j=which(tr[,"node"]==father)
        if(tr[j,"A1"] ==-1){
          tr[j,"A1"] = As[i]
        }else{
          tr[j,"A2"] = As[i]
        }
      }
    }
  }
  
  if(any(is.infinite(As))){return(list(Infinite=TRUE))}
  
  
  for(i in 2:nrow(tr)){
    father=tr[i,2]
    j=which(tr[,1]==father)
    
    if(Rs[i]==-Inf){
      A=max(1e-300,(As[j]-As[i])/As[i])
      R=-log(exp((1/eta)*log(A))+1)
      R=max(R,Rmin)
      
      Rs[i]=R
    }
    int[i]=int[j]+Rs[i]
  }
  
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,As=As,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=FALSE))
}

mcmc_eta=function(tree,epsilon,beta,ini=c(0,1),V=c(0.1,0.1),chain=NULL,
                  niter,verbose=10,silent=TRUE,Nadapt=Inf,NadaptMin=10,NadaptMax=Inf,
                  ma=-4,Ma=4,me=0.1,Me=10,proposal="bactrian",accOpt=0.3){

  Rmin=log(1e-10)-1
  maxChange=15
  ND=rank(nodes_depths_ordonnes(tree))
  tree2=build_tree(ND)
  tree$edge.length=tree2$edge.length
  
  N=tree$Nnode
  lambdaN=lambda_N(epsilon,beta,N+1)
  aN=a_N(beta,N+1)

  proba.int=function(eta,enhanced){
    p=0
    if(any(enhanced$Rs)==-Inf){return(-Inf)}
    for(i in (N+2):(2*N+1)){
      ind=which(enhanced$transform[,"node"]==i)
      son=which(enhanced$transform[,"father"]==i)
      if(length(son)>0){
        R=enhanced$Rs[son[1]]
        Rm=log(1-exp(R))
        R=max(R,Rmin)
        Rm=max(Rm,Rmin)
        Ntip=enhanced$transform[ind,"ntip"]
        ntip=enhanced$transform[son[1],"ntip"]
        cond=max(1-exp(Ntip*R)-exp(Ntip*Rm),exp(Rmin))
        p=p+(ntip)*R+(Ntip-ntip)*Rm-log(cond)
      }
    }
    return(p)
  }
  
  proba.int2=function(eta,enhanced){
    p=0
    if(any(enhanced$Rs)==-Inf){return(-Inf)}
    for(i in (N+2):(2*N+1)){
      ind=which(enhanced$transform[,"node"]==i)
      son=which(enhanced$transform[,"father"]==i)
      if(length(son)>0){
        R=enhanced$Rs[son[1]]
        Rm=log(1-exp(R))
        R=max(R,Rmin)
        Rm=max(Rm,Rmin)
        Ntip=enhanced$transform[ind,"ntip"]
        ntip=enhanced$transform[son[1],"ntip"]
        p=p+(ntip+beta)*R+(Ntip-ntip+beta)*Rm
      }else{
        son=(tree$edge[tree$edge[,1]==i,2])
        A=max(tree$tip.ab[son[1]]/tree$tip.ab[son[2]],tree$tip.ab[son[2]]/tree$tip.ab[son[1]])
        R=1/(1+A^(1/eta))
        Rm=1-R
        R=max(log(R),Rmin)
        Rm=max(log(Rm),Rmin)
        p=p+(beta+1)*R+(1+beta)*Rm
      }
    }
    return(p)
  }
  
  posterior2=function(alpha,eta,enhanced){
    if(enhanced$Infinite | eta<0){return(-Inf)}
    Y=try(aux_lik(enhanced$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    p=proba.int2(eta,enhanced)
    p2=enhanced$p_nUnsampled
    p3=enhanced$p_unsampled
    if(any(is.na(c(Y,p,p2,p3))) | any(is.nan(c(Y,p,p2,p3))) | any(is.infinite(c(Y,p,p2,p3)))){return(-Inf)}else{}
    return(Y+p-p3+p2)
  }
  
  posterior=function(alpha,eta,enhanced){
    if(enhanced$Infinite | eta<0){return(-Inf)}
    Y=try(aux_lik(enhanced$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    p=proba.int(eta,enhanced)
    if(any(is.na(c(Y,p))) | any(is.nan(c(Y,p))) | any(is.infinite(c(Y,p)))){return(-Inf)}
    return(Y+p)
  }
  
  if (is.null(chain)){
    uns=lapply(1:(N),function(i){(-1*simulate.Yi(10000,epsilon,beta,i+1))})
    a=ini[1]
    e=ini[2]
    j=1
    post_actual=-Inf
    int=enhance_eta(tree,a,beta,e,epsilon,lambdaN = lambdaN, aN=aN, uns=uns)
    post_actual=posterior(a,e,int)
    while(post_actual==-Inf){
      int=enhance_eta(tree,a,beta,e,epsilon,lambdaN = lambdaN, aN=aN, uns=uns)
      post_actual=posterior(a,e,int)
    }
    Mcmc=c(a,e,sum(int$unsampled[(N+2):(2*N+1)]),post_actual)
  }else{
    uns=chain$uns
    a=chain$a
    e=chain$e
    int=chain$int
    Mcmc=chain$mcmc
    j=nrow(chain$mcmc)
    V=chain$V
    post_actual=chain$post_actual
  }
  n=0
  for(i in 1:niter){
    if((i+j) %% Nadapt ==0 & (i+j)>NadaptMin & (i+j)<NadaptMax){ 
      if(proposal=="uniform" | proposal=="bactrian"){
        acc1=sum(diff(Mcmc[-c(1:floor(nrow(Mcmc)-Nadapt)),1])!=0)/(Nadapt-1)
        if(acc1<0.001){
          V[1]= V[1]/100
        }else if(acc1>0.999){
          V[1]= V[1]*100
        }else{
          V[1] =  V[1] * (tan((pi/2)*acc1)/tan((pi/2)*accOpt))
        }
        acc2=sum(diff(Mcmc[-c(1:floor(nrow(Mcmc)-Nadapt)),2])!=0)/(Nadapt-1)
        if(acc2<0.001){
          V[2]= V[2]/100
        }else if(acc2>0.999){
          V[2]= V[2]*100
        }else{
          V[2] =  V[2] * (tan((pi/2)*acc2)/tan((pi/2)*accOpt))
        }
      }else{
        V[1]=max(1e-10,sqrt(var(Mcmc[-c(1:floor(nrow(Mcmc)/2)),1]))*2.36)
        V[2]=max(1e-10,sqrt(var(log(Mcmc[-c(1:floor(nrow(Mcmc)/2)),2])))*2.36)
      }
      if(any(is.na(V))) V=c(1,1)
    }
    
    #Random parameter drawing : 
    if(proposal=="uniform"){
      new_a = a + (runif(1,0,1)-0.5)*3.4641016*V[1]
      new_e =exp(log(e) + (runif(1,0,1)-0.5)*3.4641016*V[2])
      
    }else if (proposal=="bactrian"){
      new_a =a + rbactrian(1)*V[1]
      new_e =exp(log(e) + rbactrian(1)*V[2])
      
    }else{
      new_a=rnorm(1,a,V[1])
      new_e=rnorm(1,e,V[2])
    }
    
    new_int=int
    post_e=posterior2(a,e,new_int)
    if(post_e>-Inf){
      ks=(N+3):(2*N+1)
      for(k in ks){
        new_new_int=change_int_eta(new_int,k,tree,a,beta,e,epsilon,lambdaN = lambdaN, aN=aN, uns=uns)
        new_post=posterior2(a,e,new_new_int)
        
        if(new_post>=post_e){
          new_int=new_new_int
          post_e=new_post
        }else{
          u=runif(1,0,1)
          if (log(u)<(new_post-post_e)){
            new_int=new_new_int
            post_e=new_post
          }
        }}
      post_new_int=posterior(a,e,new_int)
      if(post_new_int>-Inf){
        int=new_int
        post_actual=post_new_int
      }else{print("infinite")}
    }
    
    new_int=int
    if(new_a>Ma| new_a<ma){
      post=-Inf
    }else{
      post=posterior(new_a,e,int)
      
    }
    
    if(post>=post_actual){
      a=new_a
      post_actual=post
    }else{
      u=runif(1,0,1)
      if (log(u)<(post-post_actual)){
        a=new_a
        post_actual=post
      }
    }
    
    if(new_e>Me| new_e<me){
      post=-Inf
    }else{
      new_int=change_int_eta(int,NULL,tree,a,beta,new_e,epsilon,lambdaN = lambdaN, aN=aN, change=FALSE, uns=uns)
      post=posterior(a,new_e,new_int)
    }
    if(post>=post_actual){
      e=new_e
      int=new_int
      post_actual=post
    }else{
      u=runif(1,0,1)
      if (log(u)<(post-post_actual)){
        e=new_e
        int=new_int
        post_actual=post
      }
    }
    
    Mcmc=rbind(Mcmc,c(a,e,sum(int$unsampled[(N+2):(2*N+1)]),post_actual))
    
    n=n+1
    if(!silent){if(n==verbose){
      n=0
      print(paste(i,"alpha",Mcmc[i+j,1],"eta",Mcmc[i+j,2],"beta",beta,"post",Mcmc[i+j,4],"ntir",sum(int$unsampled[(N+2):(2*N+1)])))
    }}
  }
  colnames(Mcmc)=c("alpha","eta","nUnsampledled","log_post")
  return(list(mcmc=mcmc(Mcmc),post_actual=post_actual,a=a,int=int,tree=tree,V=V,e=e,uns=uns))
}


