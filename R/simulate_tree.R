lambda.epsilon=function(epsilon,beta,n,upper=10000){
  f=function(x){
    sapply(x,function(e){exp(-(beta+n+1)*e)*(1-exp(-e))^beta})
  }
  # integrate(f,epsilon,upper)
  int=2*adaptIntegrate(f, lowerLimit = epsilon, upperLimit = upper)$integral
  while(upper>0.1 & int==0){
    upper=upper/1.5
    int=2*adaptIntegrate(f, lowerLimit = epsilon, upperLimit = upper)$integral
  }
  return(int)
}

simulate.Yi=function(N,epsilon,beta,n){
  res=rep(0,N)
  if(N<1) return(c())
  
  if(beta<=0){
    g=function(x){((1-exp(-epsilon))^beta)*exp(-(beta+n+1)*epsilon)}  
    
    f=function(x){
      sapply(x,function(e){(1-exp(-e))^beta})
    }
    
    fun=function(i){
      continue=TRUE
      while(continue){
        e=rexp(1,(beta+n+1))+epsilon
        u=runif(1,0,1)
        if(u<(f(e)/g(e))){
          return(e)
          continue=FALSE
        }
      }
    }
    
    res=sapply(1:N,fun)
    
  }else{
    
    simulate.g=function(N,beta,n){
      res=rep(0,N)
      A=(beta+n+1)
      M=((beta/(A+beta))^beta)*(1+beta/A)^(-A)
      for(i in 1:N){
        u=runif(1,0,1)
        if(u<(-log(M)/(1-log(M)))){
          res[i]=runif(1,0,-log(M)/A)
        }else{
          res[i]=-log(M)/A+rexp(1,A)
        }
      }
      return(res)
    }
    
    f=function(x){
      A=beta+n+1
      M=((beta/(A+beta))^beta)*(1+beta/A)^(-A)
      sapply(x,function(e){exp(-A*e)*(1-exp(-e))^beta})
    }
    
    g=function(x){
      A=beta+n+1
      M=((beta/(A+beta))^beta)*(1+beta/A)^(-A)
      sapply(x,function(y){
        min(exp(-A*y),M)
      })
    }
    
    fun=function(i){
      continue=TRUE
      while(continue){
        e=simulate.g(1,beta,n)
        if(e>epsilon){
          u=runif(1,0,1)
          if(u<(f(e)/g(e))){
            return(e)
            continue=FALSE
          }}
      }
    }
    res=sapply(1:N,fun)
  }
  
  return(res)
}

simulate.Tau.X=function(epsilon,x,alpha,beta,n,ab=FALSE,eta=1,x.ab=1,lambda=NULL){
  if(is.null(lambda)){
    lambda=lambda.epsilon(epsilon,beta,n)
  }else{
    lambda=lambda[n]
  }
  an=sum(sapply(1:(n-1),function(i){beta(beta+1+i,beta+1+n-i)/((n+1)*beta(1+i,1+n-i))}))
  sigma=rexp(1,an)
  N=rpois(1,lambda*sigma)
  Y=c(simulate.Yi(N,epsilon,beta,n))
  Z=cumsum(c(0,Y))
  si=diff(c(0,sort(runif(N,0,sigma)),sigma))
  Tau=x^(-alpha)*(sum(exp(alpha*Z)*si))
  X=x*exp(-Z[length(Z)])
  if (ab){
    if(N==0){
      X.ab=x.ab
      return(c(Tau,X,X.ab))
    }
    
    if(eta==0){
      X.ab=x.ab/(2^N)
    }else{
      Y2=sum(log(exp(-Y*eta)+(1-exp(-Y))^eta))
      X.ab=x.ab*exp(-Z[length(Z)]*eta)/exp(Y2)
    }
    return(c(Tau,X,X.ab))
  }
  return(c(Tau,X))
}

simulate.R.K=function(beta,n){
  prob=sapply(1:(n-1),function(i){beta(beta+1+i,beta+1+n-i)/beta(1+i,1+n-i)})
  k=sample(1:(n-1),size=1,prob=prob)
  r=rbeta(1,beta+1+k,beta+1+n-k)
  return(c(r,k))
}

simulate.R=function(beta,n,K){
  r=rbeta(1,beta+1+K,beta+1+n-K)
  return(c(r,K))
}

bind.trees=function(a1,a2,root1,root2,ab=FALSE){
  if(a1$Nnode==0){return(a2)
  }else if (a2$Nnode==0){return(a1)
  }else{
    a1$root.edge=root1
    a2$root.edge=root2
    a=a1+a2
    if(ab){a$tip.ab=c(a1$tip.ab,a2$tip.ab)}
    return(a)
  }
}

simulate_tree=function(epsilon,alpha,beta,N,equal.ab=TRUE,eta=1,lambda=NULL){
  aux=function(x,n,x.ab=1){
    if(n>1){
      v1=simulate.R.K(beta,n)
      n1=v1[2]
      n2=n-n1
      x1=v1[1]*x
      x2=x-x1
      if(n1>1){
        if(equal.ab){v2.1=simulate.Tau.X(epsilon,x1,alpha,beta,n1,lambda=lambda)
        }else{v2.1=simulate.Tau.X(epsilon,x1,alpha,beta,n1,ab=TRUE,eta = eta,x.ab=(x1^eta)/(x1^eta+x2^eta),lambda=lambda)}
        tree1=aux(v2.1[2],n1)
        if(!equal.ab){
          tree1$tip.ab=v2.1[3]*tree1$tip.ab
        }
        root1=v2.1[1]
      }else{
        tree1=list(edge=matrix(c(2,1),1,2),edge.length=0,Nnode=1,tip.label=NA,tip.ab=c(1))
        if(!equal.ab){
          tree1$tip.ab=((x1^eta)/(x1^eta+x2^eta))
        }
        class(tree1)="phylo"
        
        root1=1
      }
      if(n2>1){
        if(equal.ab){v2.2=simulate.Tau.X(epsilon,x2,alpha,beta,n2,lambda=lambda)
        }else{v2.2=simulate.Tau.X(epsilon,x2,alpha,beta,n2,ab=TRUE,eta = eta,x.ab=(x2^eta)/(x1^eta+x2^eta),lambda=lambda)}
        tree2=aux(v2.2[2],n2)
        if(!equal.ab){
          tree2$tip.ab=v2.2[3]*tree2$tip.ab
        }
        root2=v2.2[1]
      }else{
        tree2=list(edge=matrix(c(2,1),1,2),edge.length=0,Nnode=1,tip.label=NA,tip.ab=c(1))
        if(!equal.ab){
          tree2$tip.ab=((x2^eta)/(x1^eta+x2^eta))
        }
        class(tree2)="phylo"
        
        root2=1
      }
      a=bind.trees(tree1,tree2,root1,root2,ab=TRUE)
      A=collapse.singles(a)
      A$tip.ab=a$tip.ab
      return(A)
    }else{
      tree=list(edge=matrix(c(2,1),1,2),edge.length=1,Nnode=1,tip.label=NA,tip.ab=1)
      class(tree)="phylo"
      return(tree)
    }
  }
  
  tree=aux(1,N)
  order=c(rep(N,N),rank(node.depth.edgelength(tree)[(N+1):(2*N-1)]))
  for(i in 1:(2*N-2)){
    tree$edge.length[i]=order[tree$edge[i,2]]-order[tree$edge[i,1]]
  }
  return(tree)
}
