\name{a_N}
\alias{a_N}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a_n
}
\description{
Values of a_n for n in 1:N
}
\usage{
a_N(beta, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
Imbalance index
}
  \item{N}{
Tip number
}
}
\value{
A vector of a_n for n in 1:N
}
\details{
\eqn{a_n = \int_(0,1) (1-r^n-(1-r)^n) r^\beta (1-r)^\beta dr}
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
