\name{simulate_yule}
\alias{simulate_yule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ranked topology with Birth-death process depths
}
\description{
Simulates the ranked topology with node depths simulated using the constant birth-death process
}
\usage{
simulate_yule(epsilon, alpha, beta, N, b, d, 
              tmax = Inf, equal.ab = TRUE, eta = 1, lambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{alpha}{
Clade age-richness index
}
  \item{beta}{
Imbalance index
}
  \item{N}{
Tip number
}
  \item{b}{
Birth rate
}
  \item{d}{
Death rate
}
  \item{tmax}{
Maximal crown age
}
  \item{equal.ab}{
If set to TRUE, all species have the same probability to go extinct first (default to TRUE)
}
  \item{eta}{
Clade abundance-richness index (if equal.ab == FALSE)
}
  \item{lambda}{
Optional, vector of lambda_epsilon returned by the function lambda_N. If set to NULL (the default) it will be computed inside the function
}
}
\value{
A phylo object with ranked shape drawn from our model, with an additional tip.ab field containing a vector of tip abundances. 
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
# Simulate a tree
set.seed(813)
tree=simulate_yule(epsilon=0.001,alpha=-1,beta=-1,
                    N=20,equal.ab=FALSE,eta=1.5, b=1, d=0.5)

# Plot the tree with dots at tips that have sizes scaling with log abundance
tree$tip.label = rep(".", length(tree$tip.label))
plot.phylo(tree, show.node.label=TRUE, 
          cex=(log(tree$tip.ab)-min(log(tree$tip.ab)-0.1))*6/diff(range(log(tree$tip.ab))), adj=0.1)
}% __ONLY ONE__ keyword per line
