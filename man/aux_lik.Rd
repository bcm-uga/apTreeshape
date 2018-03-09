\name{aux_lik}
\alias{aux_lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of the sampled node position
}
\description{
Probability of the sampled node position knowing the tree topology, the number of unsampled nodes and the interval sizes
}
\usage{
aux_lik(M, beta, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
An oject returned by the function transform
}
  \item{beta}{
Imbalance index
}
  \item{alpha}{
Clade age-richness index
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}