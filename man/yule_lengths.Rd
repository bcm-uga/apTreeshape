\name{yule_lengths}
\alias{yule_lengths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal BetaAlphaEta function
}
\description{
Simulates node depths in the birth-death model, conditionned by the number of tips and the age of the root (using the expression in "The shape and probability of reconstructed phylogenies", by Amaury Lambert and Tanja Stadler (Theoretical Population Biology, 2013))
}
\usage{
yule_lengths(N, b, d, tmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
Maximum age of the root
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}