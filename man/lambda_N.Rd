\name{lambda_N}
\alias{lambda_N}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of lambda_epsilon
}
\description{
Computes the parameter lambda_epsilon for n in 1:N
}
\usage{
lambda_N(epsilon, beta, N, upper = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
The beta parameter quantifying the balance of the tree
}
  \item{N}{
Tip number
}
  \item{upper}{
Upper bound from the integral (see lambda_n expression, default to 10000)
}
}
\details{
λ_ε = 2 ∫_(ε,∞) exp(-(β+n+1)x)(1-exp(-x))^β dx
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}


