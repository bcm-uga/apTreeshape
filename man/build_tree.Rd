\name{build_tree}
\alias{build_tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal BetaAlphaEta function
}
\description{
Builds the tree using node depths as in the birth-death model and coalescent point process;
"The shape and probability of reconstructed phylogenies", by Amaury Lambert and Tanja Stadler 
(Theoretical Population Biology, 2013)
}
\usage{
build_tree(H, tip.lab = rep(NA, length(H) + 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{
Vector of node depths
}
  \item{tip.lab}{
Vector of tip labels
}
}
\value{
A phylo object
}
\references{
Lambert, A., & Stadler, T. (2013). \emph{Birth–death models and coalescent point processes: The shape and probability of reconstructed phylogenies.} Theoretical population biology, 90, 113-128.
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}