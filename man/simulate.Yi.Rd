\name{simulate.Yi}
\alias{simulate.Yi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates the random variables Y_i
}
\description{
Simulates the random variables Y_i
}
\usage{
simulate.Yi(N, epsilon, beta, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Number of simulated Y_i
}
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
Imbalance index
}
  \item{n}{
Number of marks in the interval
}
}
\details{
\eqn{(Y_i)} are independant random variables with density \eqn{2/(\lambda_\epsilon) exp(-(\beta+n+1)x) (1-exp(-x))^\beta}
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}

