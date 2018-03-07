\name{simulate_tree}
\alias{simulate_tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates ranked topology 
}
\description{
Simulates a ranked topology with tip abundance data from our model
}
\usage{
simulate_tree(epsilon, alpha, beta, N, equal.ab = TRUE, eta = 1, lambda = NULL)
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
A phylo object with ranked shape drawn from our model. Branch lengths are so that node depths are in 1:(n-1)
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
\examples{
# Simulate a tree
set.seed(813)
tree=simulate_tree(epsilon=0.001,alpha=2,beta=-1,N=20,equal.ab=FALSE,eta=0.5)

# Plot the tree with dots at tips that have sizes scaling with log abundance
tree$tip.label = rep(".", length(tree$tip.label))
ape::plot.phylo(tree, show.node.label=TRUE, 
          cex=(log(tree$tip.ab)-min(log(tree$tip.ab)-0.1))*
          6/diff(range(log(tree$tip.ab))), adj=0.1)
}

