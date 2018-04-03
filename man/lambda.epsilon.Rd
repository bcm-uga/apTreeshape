\name{lambda.epsilon}
\alias{lambda.epsilon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes lambda_epsilon
}
\description{
Internal function needed to simulate the ranked tree shape model
}
\usage{
lambda.epsilon(epsilon, beta, n, upper = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
Imbalance index
}
  \item{n}{
Tip number
}
  \item{upper}{
Upper bound from the integral (see lambda_n expression, default to 10000)
}
}
\details{
\eqn{\lambda_\epsilon = 2 \int_(\epsilon,\infty) exp(-(\beta+n+1)x)(1-exp(-x))^\beta dx}
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
