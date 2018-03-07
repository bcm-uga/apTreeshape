\name{mcmc_eta}
\alias{mcmc_eta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference of the alpha and eta parameters
}
\description{
Run the Bayesian inference of the clade age-richness index alpha and the clade abundance richness index eta
}
\usage{
mcmc_eta(tree, epsilon, beta, ini = c(0, 1), V = c(0.1, 0.1), chain = NULL, niter, 
          verbose = 10, silent = TRUE, Nadapt = Inf, NadaptMin = 10, NadaptMax=Inf, 
          ma = -4, Ma = 4, me = 0.1, Me = 10, proposal = "bactrian", accOpt = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylo object
}
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
Imbalance index
}
  \item{ini}{
Initial alpha and eta values (default to c(0,1))
}
  \item{V}{
Initial scaling value for the mcmc proposal (default to c(0.1,0.1))
}
  \item{chain}{
Former mcmc chain (if NULL (the default), a new one is started)
}
  \item{niter}{
Number of iterations in the mcmc
}
  \item{verbose}{
Number of iterations after which the state of the mcmc is printed (if silent == FALSE)
}
  \item{silent}{
If TRUE (the default) the state of the mcmc is not printed
}
  \item{Nadapt}{
Number of iterations between each proposal scalling (default to Inf)
}
  \item{NadaptMin}{
Minimum nmber of iterations before the first proposal scalling (default to 10)
}
  \item{NadaptMax}{
Number of iterations after which the proposal stops being scalled (default to Inf)
}
  \item{ma}{
Minimal alpha value (default to -4)
}
  \item{Ma}{
Maximal alpha value (default to 4)
}
  \item{me}{
Minimal eta value (default to 0.1)
}
  \item{Me}{
Maximal eta value (default to 10)
}
  \item{proposal}{
Shape of the proposal. Can be "bactrian" (the default, ref), "uniform", or "normal"
}
  \item{accOpt}{
Optimal acceptance value (default to 0.3)
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
\examples{
seed=123
set.seed(seed)
ntip=30
tree=simulate_tree(epsilon = 0.001,alpha = -1,beta = 0,N = ntip,equal.ab = FALSE,eta =1.5)
beta=maxlik.betasplit(tree,up=10)$max_lik
extinctions = rank(tree$tip.ab)
tree$tip.label = rep(".", length(tree$tip.label))
plot.phylo(tree, show.node.label=TRUE, 
            cex=order(extinctions, seq(1,(tree$Nnode+1)))/
            ((tree$Nnode+1)/6), adj=0.1)

\dontrun{
niter=1000


chain=mcmc_eta(tree,epsilon=0.001,beta=beta,V = c(0.1,0.1),niter=niter,ini=c(0,1),
                 verbose = 50,silent = FALSE,Nadapt = 200,NadaptMin = 200)
# The initialisation of the mcmc is quiet long because 
# we begin by drawing many unsampled intervals. 
# When this is done it gets quicker. 


chain=mcmc_eta(tree,epsilon=0.001,beta=beta,niter=4000,verbose = 200,silent = FALSE,
                 chain = chain,Nadapt = 500,NadaptMin = 500,NadaptMax = 2000)

thinned=mcmc(chain$mcmc[seq(1000,5000,10),])
plot(thinned)
da=density(thinned[,"alpha"])
MPa=da$x[which.max(da$y)]
de=density(log(thinned[,"eta"]))
MPe=exp(de$x[which.max(de$y)])
print(MPa)
print(MPe)}


}% __ONLY ONE__ keyword per line
