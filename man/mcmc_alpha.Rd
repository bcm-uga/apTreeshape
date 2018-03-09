\name{mcmc_alpha}
\alias{mcmc_alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference of the alpha parameter
}
\description{
Run the Bayesian inference of the clade age-richness index alpha
}
\usage{
mcmc_alpha(tree, epsilon, beta, niter, ini = 0, V = 0.1, 
            chain = NULL, verbose = 10, silent = TRUE, Nadapt = Inf, 
            NadaptMin = 10, NadaptMax=Inf, ma = -4, Ma = 4, proposal = "bactrian", 
            accOpt = 0.3, Vmin = 0.001)
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
  \item{niter}{
Number of iterations in the mcmc
}
  \item{ini}{
Initial alpha value (default to 0)
}
  \item{V}{
Initial scaling value for the mcmc proposal (default to 0.1)
}
  \item{chain}{
Former mcmc chain (if NULL (the default), a new one is started)
}
  \item{verbose}{
Number of iterations after which the state of the mcmc is printed (if silent == FALSE)
}
  \item{silent}{
If TRUE (the default) the state of the mcmc is not printed}
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
  \item{proposal}{
Shape of the proposal. Can be "bactrian" (the default, ref), "uniform", or "normal"
}
  \item{accOpt}{
Optimal acceptance value (default to 0.3)
}
  \item{Vmin}{
Minimal scaling value for the mcmc proposal (default to 0.001)
}
}
\details{
Returns a list with a mcmc field contening the resulting chain, and other fields that are used to resume runing the inference if the chain has to be completed. 
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
\examples{
ntip=30
set.seed(123)
tree=simulate_tree(epsilon = 0.01,alpha = -1,beta = 0,N = ntip,equal.ab = TRUE)
beta=apTreeshape::maxlik.betasplit(tree,up=10)$max_lik
plot(tree)
  
niter=1000

\dontrun{
chain=mcmc_alpha(tree,epsilon=0.01,beta=beta,niter=600,V = c(0.1),ini=c(0),
                 verbose = 100,silent = FALSE,Nadapt = 100,NadaptMin = 100)
  
# Continue the same chain
chain=mcmc_alpha(tree,epsilon=0.01,beta=beta,niter=400,verbose = 100,silent = FALSE,
                 chain = chain,Nadapt = 100,NadaptMin = 500,NadaptMax = 700)
  
thinned=mcmc(chain$mcmc[seq(200,1000,10),])
plot(thinned)
da=density(thinned[,"alpha"])
MPa=da$x[which.max(da$y)]
print(MPa)

}

}% __ONLY ONE__ keyword per line
