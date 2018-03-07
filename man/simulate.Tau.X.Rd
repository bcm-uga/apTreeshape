\name{simulate.Tau.X}
\alias{simulate.Tau.X}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal BetaAlphaEta function
}
\description{
Simulates the time before the splitting of the marks into two different fragments, and the size of the fragment at this time (see appendix 1).}
\usage{
simulate.Tau.X(epsilon, x, alpha, beta, n, ab = FALSE, eta = 1, x.ab = 1, lambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{x}{
Initial fragment size
}
  \item{alpha}{
Clade age-richness index
}
  \item{beta}{
Imbalance index
}
  \item{n}{
Number of marks
}
  \item{ab}{
Boolean, should the abundances be computed as well? (default to FALSE)
}
  \item{eta}{
Clade abundance-richness index (if ab)
}
  \item{x.ab}{
Initial abundance
}
  \item{lambda}{
Optional, vector of lambda_epsilon returned by the function lambda_N. If set to NULL (the default) it will be computed inside the function
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}