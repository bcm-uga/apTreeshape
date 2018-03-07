lambda_N=function(epsilon,beta,N,upper=10000){c(0,sapply(2:N, function(x){lambda.epsilon(epsilon,beta,x,upper)}))}

a_N=function(beta,N){c(0,sapply(2:N,function(n){sum(sapply(1:(n-1),function(j){beta(beta+1+j,beta+1+n-j)/((n+1)*beta(1+j,1+n-j))}))}))}

aux_lik=function(M,beta,alpha){
  n=nrow(M)
  int=M[,"int"]
  depths=M[,"depth"]
  proba=0
  current=c(1)
  
  while (length(current)>0){
    i=which.max(depths[current])
    p=alpha*(int[current[i]])-log(sum(exp(int[current]*alpha)))
    proba=p+proba
    
    
    current=insert(current,
                   (1:n)[M[,"father"]==M[current[i],"node"]],
                   i,TRUE)
  }
  return(proba)
}

split=function(tree){
  n=tree$Nnode+1
  rep=matrix(0,3,2*n-1)
  rep[1,1:n]=rep(1,n)
  for (i in (2*n-2):1){
    rep[1,tree$edge[i,1]]=rep[1,tree$edge[i,1]]+rep[1,tree$edge[i,2]]
    if (rep[3,tree$edge[i,1]]==0){
      rep[3,tree$edge[i,1]]=rep[1,tree$edge[i,2]]
    }else{
      rep[2,tree$edge[i,1]]=rep[1,tree$edge[i,2]]
    }
  }
  return(rep)
}

insert=function(v,e,i,replace){
  if (i==0){
    c(e,v)
  }else  if(i==1){
    if (replace){
      c(e,v[-1])
    }else{c(v[1],e,v[-1])}
  }else if (i==length(v)){
    if (replace){
      c(v[-length(v)],e)
    }else{
      c(v,e)
    }
  }else if (replace){
    c(v[1:(i-1)],e,v[(i+1):length(v)])
  }else{
    c(v[1:i],e,v[(i+1):length(v)])
  }
}

rbactrian <- function(n, m=0.95){
  mBactrian = m
  sBactrian = sqrt(1-m^2)
  z = mBactrian + rnorm(n,0,1)*sBactrian
  rdunif <- runif(n,0,1)<0.5
  sign <- ifelse(rdunif,-1,1)
  z=z*sign
  return(z)
}

transform=function(tree,unsampled,list_depth=NULL,change_depth=0){
  split=split(tree)
  ntip=tree$Nnode+1
  depth=c(rep(0,ntip),depth(tree))
  p_unsampled=0
  compute_depth=is.null(list_depth)
  if(compute_depth){list_depth=list(c())}
  fils=tree$edge[tree$edge[,1]==(ntip+1),2]
  M=c(ntip+1,-1,depth[ntip+1],ntip)
  k=2*ntip
  id.unsampled=1
  for (i in (ntip+2):(2*ntip-1)){
    father=tree$edge[tree$edge[,2]==i,1]
    node_father=father
    pf=depth[i]
    pp=depth[father]
    fils=tree$edge[tree$edge[,1]==i,2]
    if (unsampled[i]>0){
      if(compute_depth | (change_depth==i)){
        depths=runif(unsampled[i],pf,pp)
        depths=depths[order(depths,decreasing = TRUE)]
        list_depth[[i]]=depths
      }else{
        depths=list_depth[[i]]
      }
      for (j in 1:unsampled[i]){
        M=rbind(M,c(k,father,depths[j],split[1,i]))
        father=k
        k=k+1
        id.unsampled=id.unsampled+1
      }
      p_unsampled=p_unsampled+sum(log(1:unsampled[i]))-log(pp-pf)*unsampled[i]
    }
    M=rbind(M,c(i,father,pf,split[1,i]))
  }
  colnames(M)=c("node","father","depth","ntip")
  for(j in 1:M[1,"depth"]){
    N=sum(M[,"depth"]>(j-1) & M[,"depth"]<(j))
    if(N>0){
      p_unsampled=p_unsampled-sum(log(1:N))
    }
  }
  
  return(list(M=M,p_unsampled=p_unsampled,list_depth=list_depth))
}

enhance=function(tree,alpha,beta,epsilon,lambdaN,aN){
  N=tree$Nnode
  unsampled=c(-1,rep(-1,2*N))
  p_nUnsampled=0
  
  for (i in (N+3):(2*N+1)){
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
  
  tr=transform(tree,unsampled)
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  
  int=c(0,rep(-Inf,nrow(tr)-1))
  Rs=c(1,rep(-1,nrow(tr)-1))
  
  for(i in 2:nrow(tr)){
    if(int[i]==-Inf){
      father=tr[i,2]
      j=(1:nrow(tr))[tr[,1]==father]
      if(father<(2*N+2)){
        R=simulate.R(beta,tr[j,4],tr[i,4])[1]
        int[i]=int[j]+log(R)
        Rs[i]=R
        if(i<nrow(tr)){
          k=((i+1):nrow(tr))[tr[(i+1):nrow(tr),2]==father]
          int[k]=int[j]+log(1-R)
          Rs[k]=(1-R)
        }
      }else{
        R=exp(-1*simulate.Yi(1,epsilon,beta,tr[i,4]))
        int[i]=int[j]+log(R)
        Rs[i]=R
      }
    }
  }
  
  
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=FALSE))
}

change_int=function(enhanced,ind,tree,alpha,beta,epsilon,lambdaN,aN){
  changeMax=5
  N=tree$Nnode
  unsampled=enhanced$unsampled
  int=enhanced$int
  p_nUnsampled=enhanced$p_nUnsampled
  list_depth=enhanced$list_depth
  continue=TRUE
  if(is.null(ind) ){
    rm=c()
  }else if (ind==(N+2)) {
    rm=c()
  }else{
    n=extract.clade(tree,ind)$Nnode+1
    rm=c()
    tr=enhanced$transform
    node=which(tr[,1]==ind)
    if(unsampled[ind]>0){
      rm=(node-1):(node-unsampled[ind])
      node=node-unsampled[ind]
    }
    
    if(n==1){unsampled[ind]=-1
    }else if(ind==(N+2)){unsampled[ind]=-1}else{
      lambda=lambdaN[n]
      an=aN[n]
      p_nUnsampled=p_nUnsampled+(unsampled[ind]+1)*log(lambda+an)-(unsampled[ind])*log(lambda)
      u=runif(1)
      if(u<(0.4)){
        unsampled[ind]=unsampled[ind]-1
      }else if(u>0.6){
        unsampled[ind]=unsampled[ind]+1
      }
      unsampled[ind]=unsampled[ind]+floor(rnorm(1,mean = 0.5,sd = 10))
      if(unsampled[ind]<0){
        proba.unsampled=-Inf
        continue=FALSE
      }else{
        p_nUnsampled=p_nUnsampled-(unsampled[ind]+1)*log(lambda+an)+(unsampled[ind])*log(lambda)
        continue=TRUE
      }
      if(is.na(unsampled[ind])) print(paste("a",alpha,"b",beta,"t",unsampled[ind],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  
  if(! continue){
    enhanced$p_nUnsampled=-Inf
    return(enhanced)
  }
  
  tr=transform(tree,unsampled,list_depth = list_depth,change_depth = max(ind,0))
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  oldRs=c()
  
  if(length(rm)>0){
    oldRs=enhanced$Rs[rm]
    Rs=enhanced$Rs[-rm]
    depth=enhanced$transform[-rm,"depth"]
    
  }else{
    
    Rs=enhanced$Rs
    depth=enhanced$transform[,"depth"]
    
  }
  if(is.null(ind)){
    Rs=enhanced$Rs
  }else{if(unsampled[ind]>-1){
    add=c(-1,sample(c(oldRs,rep(-1,changeMax+max(0,unsampled[ind]-length(oldRs)))),unsampled[ind]))
    if(node<length(Rs)){
      Rs=c(Rs[1:(node-1)],add,Rs[(node+1):length(Rs)])
      depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"],depth[(node+1):length(depth)])
      
    }else{
      Rs=c(Rs[1:(node-1)],add)
      depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"])}
  }}
  
  tr[,"depth"]=depth
  
  int=c(0,rep(-Inf,nrow(tr)-1))
  if(! is.null(ind)){if(ind>(N+2)){father=tr[node,2]
  Ks=which(tr[,2]==father)
  j=(1:nrow(tr))[tr[,1]==father]
  if(Rs[Ks[1]]==1){R=simulate.R(beta,tr[j,4],tr[Ks[1],4])[1]
  Rs[Ks[1]]=R
  if(length(Ks)==2){Rs[Ks[2]]=(1-R)}}}}
  Rs[1]=1
  
  for(i in 2:nrow(tr)){
    if(int[i]==-Inf){
      father=tr[i,2]
      j=(1:nrow(tr))[tr[,1]==father]
      if(father<(2*N+2)){
        if(Rs[i]==-1){
          R=simulate.R(beta,tr[j,4],tr[i,4])[1]
          Rs[i]=R}
        R=Rs[i]
        int[i]=int[j]+log(R)
        if(i<nrow(tr)){
          k=((i+1):nrow(tr))[tr[(i+1):nrow(tr),2]==father]
          int[k]=int[j]+log(1-R)
          Rs[k]=(1-R)
        }
      }else{
        if(Rs[i]==-1){
          R=exp(-1*simulate.Yi(1,epsilon,beta,tr[i,4]))
          Rs[i]=R}
        R=Rs[i]
        int[i]=int[j]+log(R)
        
      }
    }
  }
  if(any(Rs<0) ) {print(ind); debug(change_int); change_int(enhance,ind,tree,alpha,beta,epsilon)}
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=FALSE))
}

mcmc_alpha=function(tree,epsilon,beta,niter,ini=0,V=0.1,chain=NULL,
                    verbose=10,silent=TRUE,Nadapt=Inf,NadaptMin=10,NadaptMax=Inf,
                    ma=-4,Ma=4,proposal="bactrian",accOpt=0.3,Vmin=0.001){
  
  
  ND=rank(nodes_depths_ordonnes(tree))
  tree2=build_tree(ND)
  tree$edge.length=tree2$edge.length
  N=tree$Nnode
  lambdaN=lambda_N(epsilon,beta,N+1)
  aN=a_N(epsilon,beta,N+1)
  ntip=tree$Nnode+1
  ordre=rank(max(depth(tree))+1-depth(tree))
  X=rep(1,ntip-1)
  for (i in 1:(ntip-1)){X[ordre[i]]=((ntip+1):(ntip*2-1))[i]}
  
  
  posterior2=function(alpha,enhance){
    if(enhance$Infinite){return(-Inf)}
    Y=try(aux_lik(enhance$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    p2=enhance$p_nUnsampled
    p3=enhance$p_unsampled
    if(any(is.na(c(Y,p2,p3))) | any(is.nan(c(Y,p2,p3))) | any(is.infinite(c(Y,p2,p3)))){return(-Inf)}
    return(Y-p3+p2)
  }
  
  posterior=function(alpha,enhance){
    if(enhance$Infinite ){return(-Inf)}
    Y=try(aux_lik(enhance$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    if(is.na(Y) | is.nan(Y) | is.infinite(Y)){return(-Inf)}
    return(Y)
  }
  
  if (is.null(chain)){
    a=ini
    j=1
    post_actuel=-Inf
    while(post_actuel==-Inf){
      int=enhance(tree,a,beta,epsilon,lambdaN,aN)
      post_actuel=posterior(a,int)
      post_actuel_int=posterior2(a,int)
    }
    Mcmc=c(a,sum(int$unsampled[-(1:(N+2))]),post_actuel)
  }else{
    a=chain$a
    int=chain$int
    Mcmc=chain$mcmc
    j=nrow(chain$mcmc)
    V=chain$V
    post_actuel=chain$post_actuel
    post_actuel_int=posterior2(a,int)
  }
  n=0
  for(i in 1:niter){
    if((i+j) %% Nadapt ==0 & (i+j)>NadaptMin & (i+j)<NadaptMax){ 
      if(proposal=="uniform" | proposal=="bactrian"){
        acc=sum(diff(Mcmc[-c(1:floor(nrow(Mcmc)-Nadapt)),1])>0)/(Nadapt-1)
        if(acc<0.001){
          V=V/100
        }else if(acc>0.999){
          V=V*100
        }else{
          V = V * (tan((pi/2)*acc)/tan((pi/2)*accOpt))
        }
        V=max(V,Vmin)
      }else{
        V=max(Vmin,sqrt(var(Mcmc[-c(1:floor(nrow(Mcmc)/2)),1]))*2.36)
      }
      if(is.na(V)) V=1
    }
    
    #On tire les paramètres proposés : 
    post_actuel=try(posterior(a,int))
    if(proposal=="uniform"){
      new_a = a + (runif(1,0,1)-0.5)*3.4641016*V
    }else if (proposal=="bactrian"){
      new_a =a + rbactrian(1)*V
    }else{
      new_a=rnorm(1,a,V)}
    if(new_a>Ma| new_a<ma){
      post=-Inf
    }else{
      post=try(posterior(new_a,int))}
    
    if(post>=post_actuel){
      a=new_a
      post_actuel=post
      post_actuel_int=posterior2(a,int)
    }else{
      u=runif(1,0,1)
      if (log(u)<(post-post_actuel)){
        a=new_a
        post_actuel=post
        post_actuel_int=posterior2(a,int)
      }
    }
    
    ks=(N+3):(2*N+1)
    for(k in ks){
      new_int=change_int(int,k,tree,a,beta,epsilon,lambdaN,aN)
      post=posterior2(a,new_int)
      
      if(post>=post_actuel_int){
        int=new_int
        post_actuel_int=post
        post_actuel=posterior(a,int)
      }else{
        u=runif(1,0,1)
        if (log(u)<(post-post_actuel_int)){
          int=new_int
          post_actuel_int=post
          post_actuel=posterior(a,int)
        }
      }}
    
    
    Mcmc=rbind(Mcmc,c(a,sum(int$unsampled[-(1:(N+2))]),post_actuel))
    
    n=n+1
    if(n==verbose){
      n=0
      print(paste(i,"alpha",Mcmc[i+j,1],"beta",beta,"post",Mcmc[i+j,3]))
    }
  }
  colnames(Mcmc)=c("alpha","nUnsampled","log_post")
  return(list(mcmc=mcmc(Mcmc),post_actuel=post_actuel,a=a,int=int,tree=tree,V=V))
}

